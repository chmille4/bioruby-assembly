module Bio
  class Assembly
    class Maf < Bio::Assembly
      
      # register parser with superclass
      register_parser :maf

      def initialize(path)
         @file = File.new(path, 'r') 
         # TO DO
      end

    end # end Maf
  end # end Assembly  
end # end Bio