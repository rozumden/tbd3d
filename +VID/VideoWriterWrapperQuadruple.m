classdef VideoWriterWrapperQuadruple < handle
   properties
      fps = 5
      repeat = 1
      path
      video
      first
      img_buffer = []
      img_buffer2 = []
      img_buffer3 = []
   end

   methods
      function this = VideoWriterWrapperQuadruple(folder, name, varargin)
         this = cmp_argparse(this, varargin);
         [~,b,~] = fileparts(name);
         this.path = fullfile(folder,'results',[b '_result.avi']);
         this.video = VideoWriter(this.path);
         this.video.FrameRate = this.fps;
         open(this.video);
         this.first = true;
      end


      function [] = write_buffer_img(this, img)
         this.img_buffer = img;
      end

      function [] = write_buffer2_img(this, img)
         this.img_buffer2 = img;
      end    

      function [] = write_buffer3_img(this, img)
         this.img_buffer3 = img;
      end        

      function [] = write_img(this, img)
         this.img_buffer = imresize(this.img_buffer, [size(img,1) NaN]);
         this.img_buffer2 = imresize(this.img_buffer2, [size(img,1) NaN]);
         this.img_buffer3 = imresize(this.img_buffer3, [size(img,1) NaN]);
         img = cat(2, this.img_buffer, this.img_buffer2, this.img_buffer3, img);
         
         fctr = 1;
        
         for tt = 1:fctr*this.repeat
            writeVideo(this.video,uint8(img));
         end
      end

      function [] = close(this)
         close(this.video);
      end
   end

end