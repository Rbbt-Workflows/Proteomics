module Proteomics
  FOLD_SEP = "··"
  FOLD_SEP2 = "·##·"
  def self.unfold_traverse(stream, workflow, task, input = nil, other_inputs = {}, &block)
    key_field, unfold_field = Misc.process_options other_inputs, :key_field, :unfold_field
    unfolded_stream = Misc.open_pipe do |sin|
      TSV.traverse stream, :into => sin do |key,values|
        key = key.first if Array === key
        num = 0
        total = values.length.to_s
        res = values.collect{|v| num += 1; v + FOLD_SEP + num.to_s + FOLD_SEP + total + FOLD_SEP + key.gsub(FOLD_SEP, FOLD_SEP2) }
        res.extend MultipleResult
        res
      end.join
      sin.close
    end

    job = workflow.job(task, nil, other_inputs.merge(input => unfolded_stream))

    parser = TSV::Parser.new job.exec(:stream)
    options = parser.options
    options[:key_field] = key_field if key_field
    options[:fields] = [unfold_field] + options[:fields] if unfold_field
    dumper = TSV::Dumper.new options
    dumper.init

    current_values = {}
    last_pos = parser.options[:fields].length
    TSV.traverse parser, :into => dumper  do |key,values|
      values[last_pos-1] ||= []
      key = key.first if Array === key

      unfold_key, num, total, folded_key = key.split(FOLD_SEP)

      folded_key = folded_key.gsub(FOLD_SEP2, FOLD_SEP)

      values = [unfold_key] + values if unfold_field

      values = yield folded_key, unfold_key, values if block_given?

      values = values.collect do |v|
        case v
        when nil
          nil
        when Array
          v.collect{|e| e.to_s.gsub(';', ',').gsub('|', ';') }
        when String
          v.gsub(';', ',').gsub('|', ';') 
        else
          v
        end
      end

      current_values[folded_key] ||= []
      current_values[folded_key][num.to_i - 1] = values

      next unless current_values[folded_key].compact.length == total.to_i
        
      res = [folded_key, Misc.zip_fields(current_values[folded_key])]

      current_values.delete(folded_key)

      res
    end

    dumper
  end

end
