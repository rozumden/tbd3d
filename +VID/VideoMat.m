classdef VideoMat < VID.IVideo
   properties(Access = private)
      k
      sz 
      w
      h
   end

   properties(Access = public)
      mat
   end

   methods
      function this = VideoMat(mat)
         this.mat = mat;
         this.k = 1;
         this.sz = size(mat, 4);
         [this.h,this.w,~] = size(mat);
      end

      function [] = play(this, PSF) 
         VkPOS = this.mat;
         VkPOS(:,:,1 , : ) = squeeze(VkPOS(:,:,1,:))+255*uint8(PSF>0);
         implay(VkPOS);
      end 

      function [matout] = get_mat(this)
         matout = this.mat;
      end

      function frame = get_next(this)
         frame = [];
         if this.has_next()
            frame = this.get_frame(this.k);
            this.k = this.k + 1;
         end
      end

      function has = has_next(this)
         has = this.k <= this.sz;
      end

      function [framed, frame] = get_frame(this, n)
         frame = this.mat(:,:,:,n);
         framed = im2double(frame);
      end

      function num = length(this)
         num = this.sz;
      end

      function num = width(this)
         num = this.w;
      end

      function num = height(this)
         num = this.h;
      end
      
      function [] = close(this)

      end

   end

end

