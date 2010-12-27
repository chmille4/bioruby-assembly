module Bio
  class Assembly
    class Read
      
      module Ace 
        attr_accessor :base_sequences
        
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