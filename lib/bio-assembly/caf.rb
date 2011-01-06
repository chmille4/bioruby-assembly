module Bio
  class Assembly
    class Caf < Bio::Assembly
      
      # register parser with superclass
      register_parser :caf

      def initialize(path)
         @file = File.new(path, 'r') 
         
      end

    end # end Caf
  end # end Assembly  
end # end Bio