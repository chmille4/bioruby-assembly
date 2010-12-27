module Bio
   class Assembly
      class Contig
      
         # ace specific methods for contig objects
         module Ace
         
            def to_ace
              ace = ""
              ace += ['CO', name, num_bases, num_reads, num_base_segments, orientation].join(' ') + "\n"
              ace += seq.to_s.gsub(Regexp.new(".{1,50}"), "\\0\n") + "\n"
              ace += "BQ\n"
              last_stop = quality.size - 1
              (quality.size/50+1).times do |i|
                start = i * 50
                stop = (i+1) * 50 - 1
                stop = last_stop if stop > last_stop
                ace += ' ' + quality[start..stop].join(' ')  + "\n"
              end
              ace += "\n"

              # holds BS data for reads
              bs_str = ""
              # holds RD, QA, and DS data for reads 
              rest_str = ""
              @reads.values.sort.each do |read|
                ace += read.to_ace_af
                bs_str += read.to_ace_bs
                rest_str += read.to_ace_rest
              end

              # compile data in correct order
              ace += bs_str
              ace += "\n"
              ace += rest_str
              ace
            end
            
         end
         
      end
   end
end