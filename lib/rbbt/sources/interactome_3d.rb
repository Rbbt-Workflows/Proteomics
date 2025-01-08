require 'rbbt-util'

module Interactome3d
  extend Resource
  self.subdir = 'share/databases/interactome3d'

  self.claim self.proteins_tsv, :proc do |filename|
    tsv = TSV.open('http://interactome3d.irbbarcelona.org/user_data/human/download/complete/proteins.dat', 
                   :header_hash => '',
                   :merge => true)
    tsv.to_s
  end

  self.claim self.interactions_tsv, :proc do |filename|
    tsv = TSV.open('http://interactome3d.irbbarcelona.org/user_data/human/download/complete/interactions.dat', 
                   :header_hash => '',
                   :merge => true)
    tsv.to_s
  end

  self.claim self.proteins, :proc do |filename|
    self.proteins_tsv.produce
    Open.link self.proteins_tsv.find, filename
  end
  
  self.claim self.interactions, :proc do |filename|
    self.interactions_tsv.produce
    Open.link self.interactions_tsv.find, filename
  end


  self.claim self[".source"], :rake do
    rule /(proteins|interactions)_\d+\.tgz/ do |t|
      url = "https://interactome3d.irbbarcelona.org/user_data/human/download/complete/#{File.basename(t.name)}"
      TmpFile.with_file do |tmpfile|
        Open.download(url, tmpfile)
        Open.mv(tmpfile, t.name)
      end
    end
  end

  self.claim self.pdbs, :proc do |directory|
    list = Open.read("https://interactome3d.irbbarcelona.org/downloadset.php?queryid=human&release=current&path=complete")
    proteins = list.scan(/proteins_\d+\.tgz/).uniq
    interactions = list.scan(/interactions_\d+\.tgz/).uniq

    proteins.each do |file|
      Misc.untar(Interactome3d[".source"][file].produce, directory)
    end

    interactions.each do |file|
      Misc.untar(Interactome3d[".source"][file].produce, directory)
    end
  end

  #self.claim self["proteins_bulk.tar.gz"], :proc do |filename|
  #  list = Open.read("https://interactome3d.irbbarcelona.org/downloadset.php?queryid=human&release=current&path=complete")
  #  proteins = list.scan(/proteins_\d+\.tgz/).uniq

  #  TmpFile.with_dir do |tmpdir|
  #    proteins.each do |file|
  #      Misc.untar(Interactome3d[".source"][file].produce, tmpdir)
  #    end
  #    Misc.tarize(tmpdir, filename)
  #  end
  #end

  #self.claim self["interactions_bulk.tar.gz"], :proc do |filename|
  #  list = Open.read("https://interactome3d.irbbarcelona.org/downloadset.php?queryid=human&release=current&path=complete")
  #  interactions = list.scan(/interactions_\d+\.tgz/).uniq

  #  TmpFile.with_dir do |tmpdir|
  #    interactions.each do |file|
  #      Misc.untar(Interactome3d[".source"][file].produce, tmpdir)
  #    end
  #    Misc.tarize(tmpdir, filename)
  #  end
  #end

  def self.pdb(pdb)
    raise "Could not produce Interactome3d pdbs" unless Interactome3d.pdbs.produce
    Interactome3d.pdbs[pdb].find
  end
end

if __FILE__ == $0
  Log.severity = 0
  begin
    Interactome3d["proteins_bulk.tar.gz"].produce
    Interactome3d["interactions_bulk.tar.gz"].produce
  rescue
    sleep 5
    retry
  end
end
