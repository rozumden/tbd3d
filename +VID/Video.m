classdef Video < VID.IVideo
   properties(Access = private)
      file
      object
      k
   end

   methods
      function this = Video(folder, name, datasetname)
         if ~exist('datasetname')
            datasetname = 'seq';
         end
         file = fullfile(folder,datasetname,name);
         this.file = file;
         if exist(file) == 7 || exist(file) == 0 
            if ~any(exist(name) == [0 7])
               file = name;
            else
               error('No such file!');
               return
            end
         end

         global check_file;
         if isempty(check_file)
            check_file = containers.Map;
         end
         if check_file.isKey(file)
            v = check_file(file);
            if ~isvalid(v)
               v = VideoReader(file);
               check_file(file) = v;
            end
         else
            was_caught = false;
            try
               v = VideoReader(file);
            catch 
               v = VideoReader(file);
            end
            check_file = [check_file; containers.Map(file, v)];
         end
         this.object = v;
         this.k = 1;
      end

      function [] = close(this)
         this.object.close();
      end

      function frame = get_next(this)
         frame = [];
         if this.has_next()
            frame = this.get_frame(this.k);
            this.k = this.k + 1;
         end
      end

      function has = has_next(this)
         has = this.k <= this.length();
      end

      function [framed, frame] = get_frame(this, n)
         frame = this.object.read(n);
         framed = im2double(frame);
      end

      function num = length(this)
         num = this.object.NumberOfFrames;
      end

      function num = width(this)
         num = this.object.Width;
      end

      function num = height(this)
         num = this.object.Height;
      end
      
   end

end

