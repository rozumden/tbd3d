classdef VideoWriterWrapperDouble < handle
   properties
      fps = 5
      repeat = 1
      path
      video
      first
      img_buffer = []
   end

   methods
      function this = VideoWriterWrapperDouble(folder, name, varargin)
         this = cmp_argparse(this, varargin);
         [~,b,~] = fileparts(name);
         this.path = fullfile(folder,'results',[b '_result.avi']);
         this.video = VideoWriter(this.path);
         this.video.FrameRate = this.fps;
         open(this.video);
         this.first = true;
      end

      function img = get_img(this)
         set(gca,'position',[0 0 1 1],'units','normalized');
         set(gca,'xtick',[])
         set(gca,'ytick',[])
         set(gca, 'box','off')
         set(gca,'Visible','off')
         axis equal
         F = getframe(gcf);
         [img,~] = frame2im(F);
         w = size(img,2);
         w1 = ceil(w/10);
         img = img(:,30:(w-w1),:);
      end

      function [] = write_buffer(this)
         this.img_buffer = this.get_img();
      end   

      function [] = write(this)
         img = this.get_img();
         img = cat(2, this.img_buffer, img);
         
         fctr = 1;
         % if this.first 
         %    this.first = false;
         %    fctr = this.fps;
         % end

         for tt = 1:fctr*this.repeat
            writeVideo(this.video,uint8(img));
         end
      end

      function [] = close(this)
         img = this.get_img();
         % for kk = 1:this.fps*this.repeat
            % writeVideo(this.video,uint8(img));
         % end
         close(this.video);
      end
   end

end