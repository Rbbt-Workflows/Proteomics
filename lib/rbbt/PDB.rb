require 'rbbt-util'
require 'rbbt/PDB/distances'

module PDB
  USE_REPOSITORY_DIR = true

  def self.cache(code)
    Rbbt.var.Proteomics.PDB[code.to_s]
  end

  def self.pdb_stream(pdb = nil, pdbfile = nil)
    return StringIO.new(pdbfile) if (pdb.nil? or pdb.empty?) and not pdbfile.nil? and not pdbfile.empty?
    raise "No PDB specified" if pdb.nil? or pdb.empty?
    pdb = pdb.strip

    if Open.exists?(pdb)
      Open.open(pdb)
    else
      begin
        content = Persist.persist("PDB file", :text, :persist => true, :dir => cache(:pdb_files), :other => {:pdb => pdb}) do 
          Misc.insist(1) do
            if pdb.length > 5 || Open.remote?(pdb)
              Open.read(pdb)
            else
              pdb = pdb.sub(/^=/,'')
              Open.read("http://www.rcsb.org/pdb/download/downloadFile.do?fileFormat=pdb&compression=NO&structureId=#{pdb}")
            end
          end
        end
        StringIO.new(content)
      rescue Exception
        raise "No valid pdb provided: #{ Misc.fingerprint pdb }"
      end
    end
  end

  def self.pdb_chain_sequences(pdb = nil, pdbfile = nil)
    Persist.persist("Chain sequences", :marshal, :persist => true, :dir => cache(:chain_sequences), :other => {:pdb => pdb, :pdbfile => pdbfile}) do 
      
      chains = {}
      stream = pdb_stream(pdb, pdbfile)
      TSV.traverse stream, :type => :array do |line|
        break if line =~ /^END/
        next unless line =~ /^ATOM/ and line.length > 21
        chain = line[20..21].strip
        aapos = line[22..25].to_i
        aa    = line[17..19]

        next if aapos <= 0

        chain = "A" if chain == ""
        chains[chain] ||= Array.new
        chains[chain][aapos-1] = aa
      end
      stream.close unless stream.closed?

      chains.each do |chain,chars|
        chains[chain] = chars.collect{|aa| aa.nil? ? '?' : Misc::THREE_TO_ONE_AA_CODE[aa.downcase] || aa.strip} * ""
      end

      chains
    end
  end

  def self.pdb_alignment_map(sequence, pdb = nil, pdbfile = nil)
    chains = pdb_chain_sequences(pdb, pdbfile)
    Persist.persist("PDB aligment map", :marshal, :persist => true, :dir => cache(:pdb_alignment_map), :other => {:pdb => pdb, :pdbfile => pdbfile, :sequence => sequence}) do 
      map = {} 
      chains.each do |chain,chain_sequence|
        map[chain] = SmithWaterman.alignment_map(sequence, chain_sequence)
      end
      map
    end
  end

end
