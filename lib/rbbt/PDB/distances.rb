module PDB

  MAX_DISTANCE = 10

  def self.atom_distance(pdb = nil, pdbfile = nil, max_distance = MAX_DISTANCE)
    Persist.persist("Atom distances", :marshal, :dir => cache(:atom_distance), :other => {:pdb => pdb, :pdbfile => pdbfile, :max_distance => max_distance}) do 
      Log.low "Computing atom distances (#{ max_distance }): #{pdb || "pdbfile"}"
      stream = pdb_stream(pdb, pdbfile)

      atom_positions = {}
      while line = stream.gets
        break if line =~ /^END/
        #next unless line =~ /^ATOM/

        # THIS LINE FIXES A PROBLEM WITH INTERACTOME 3D PDBs
        # THAT BREAK THE FIRST LINE AFTER A 'TER'
        line = line.sub(/^\d+/,'')

        next unless line =~ /^\w*ATO?M/ and line.length > 53
        code = line[13..26]
        x = line[30..37].to_f
        y = line[38..45].to_f
        z = line[46..53].to_f
        num = code[9..13].to_i

        atom_positions[code] = [x,y,z,num]
      end

      atom_positions

      atoms = atom_positions.
        sort_by{|a,values| values[3] }.
        collect{|atom,v| atom }

      atom_distances = []

      while atom1 = atoms.shift
        position1 = atom_positions[atom1]

        atoms.each do |atom2|
          position2 = atom_positions[atom2]
          next if (position1[3] == position2[3]) 
          next if ((position1[3] == position2[3] + 1) and 
                   ((atom1[0] == "C" and atom2[0] == "N") or 
                    (atom1[0] == "0" and atom2[0] == "N") or 
                    (atom1[0] == "C" and atom2[0..1] == "CA")))

          position2 = atom_positions[atom2]
          dx = position1[0] - position2[0]
          dy = position1[1] - position2[1]
          dz = position1[2] - position2[2]

          next if dx.abs > max_distance or dy.abs > max_distance or dz.abs > max_distance
          dist = Math.sqrt(dx**2 + dy**2 + dz**2)
          next if dist > max_distance
          atom_distances << [atom1, atom2, dist]
        end

      end

      atom_distances
    end
  end

  def self.close_residues(distance, pdb = nil, pdbfile = nil)
    atom_distances = atom_distance(pdb, pdbfile, [distance, MAX_DISTANCE].max)

    close_residues = {}
    Log.low "Computing residue distances (#{distance}): #{pdb || "pdbfile"}"
    atom_distances.each do |atom1, atom2, dist|
      next if dist > distance
      aa1 = atom1.strip.split(/\s+/).values_at(2,3) * ":"
      aa2 = atom2.strip.split(/\s+/).values_at(2,3) * ":"
      aa1 = "A" + aa1 if aa1 =~ /^\d/
      aa2 = "A" + aa2 if aa2 =~ /^\d/
      aa1 = aa1[0] + ":" + aa1[1..-2] if aa1[-1] == ":"
      aa2 = aa2[0] + ":" + aa2[1..-2] if aa2[-1] == ":"
      close_residues[aa1] ||= []
      close_residues[aa1] << aa2
    end

    close_residues.each do |aa1, list| list.uniq! end

    close_residues
  end

  def self.interface_residues(pdb = nil, pdbfile = nil, distance = 8)
    Persist.persist("Interface residues", :marshal, :dir => cache(:interface_residues), :other => {:pdb => pdb, :pdbfile => pdbfile, :distance => distance}) do 
      close_residues = PDB.close_residues(distance, pdb, pdbfile)

      interfaces = []
      close_residues.each do |r1,other|
        c1, p1 = r1.split(":")
        other.each do |r2|
          c2, p2 = r2.split(":")
          next if c1 == c2
          interfaces << [r1, r2]
        end
      end

      interfaces
    end
  end
end
