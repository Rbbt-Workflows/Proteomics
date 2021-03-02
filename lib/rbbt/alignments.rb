require 'rbbt/tools/ssw'

module Proteomics
  def self.alignment_map(source, target)
    alignment_source, alignment_target = SmithWaterman.align(source, target)
    map = {}

    offset_source, alignment_source = alignment_source.match(/^(_*)(.*)/).values_at( 1, 2)
    offset_target, alignment_target = alignment_target.match(/^(_*)(.*)/).values_at( 1, 2)

    gaps_source = 0 
    gaps_target = 0
    miss_match = 0
    alignment_source.chars.zip(alignment_target.chars).each_with_index do |p,i|
      char_source, char_target = p
      gaps_source += 1 if char_source == '-'
      gaps_target += 1 if char_target == '-'
      source_pos = i + 1 + offset_source.length - gaps_source
      target_pos = i + 1 + offset_target.length - gaps_target
      if char_source != char_target or char_source == "-"
        miss_match += 1
      else
        map[source_pos] = target_pos 
      end
    end

    if miss_match + gaps_source > alignment_source.length.to_f / 2
      {}
    else
      map
    end
  end
end
