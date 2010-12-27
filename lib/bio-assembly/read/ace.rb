module Bio
  class Assembly
    class Read
      
      # ace specific methods for read objects
      module Ace 
        attr_accessor :base_sequences
        
        def to_ace
          ace += ""
          # holds BS data for reads
          bs_str = ""
          # holds RD, QA, and DS data for reads 
          rest_str = ""
          ace += to_ace_af
          bs_str += to_ace_bs
          rest_str = to_ace_rest

          # compile data in correct order
          ace += bs_str
          ace += "\n"
          ace += rest_str
          ace
        end

        def to_ace_bs
          bs_str = ""
          unless base_sequences.nil?
            base_sequences.each do |bs|
              bs_str += ['BS', bs.from, bs.to, bs.read_name].join(' ') + "\n"
            end
          end
          bs_str
        end

        def to_ace_af
          ['AF', name, orientation, from].join(' ') + "\n"
        end

        def to_ace_rest
          rest_str = ""
          rest_str += ['RD', name, num_bases, 0, 0].join(' ') + "\n"
          rest_str += seq.to_s.gsub(Regexp.new(".{1,50}"), "\\0\n")  + "\n"
          rest_str += ['QA', clear_range_from, clear_range_to, clear_range_from, clear_range_to].join(' ') + "\n"
          rest_str += ['DS', 'CHROMAT_FILE:', name, 'PHD_FILE:', "#{name}.phd.1", 'TIME:', Time.now].join(' ') + "\n"
          rest_str
        end
        
        def add_base_sequence(from, to, read_name)
          @base_sequences = Array.new if @base_sequences.nil?
          @base_sequences.push BaseSequence.new(from, to, read_name)
        end
        
        class BaseSequence
          attr_accessor :from, :to, :read_name
          
          def initialize(from, to, read_name)
            @from = from
            @to = to
            @read_name = read_name
          end
          
          def <=>(other)
            unless other.kind_of?(Bio::Assembly::Read::BaseSequence)
              raise "[Error] markers are not comparable"
            end
            if self.from == other.from
              # sort by to if froms are identical
              return self.to.<=>(other.to)
            else
              return self.from.<=>(other.from)
            end
          end
          
        end
        
      end
    end
  end
end