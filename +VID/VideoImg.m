classdef VideoImg < VID.IVideo
   properties
      digits = 8
      imgfolder
      k
      sz 
      w
      h
      mask
      apply_function = []
   end

   methods
      function this = VideoImg(folder, name, datasetname, varargin)
         this = cmp_argparse(this, varargin);
         if ~exist('datasetname')
            datasetname = 'seq_img';
         end
         [~,nameseq,~] = fileparts(name);
         file = fullfile(folder,datasetname,nameseq);
         this.imgfolder = file;
         if exist(file) == 0 
            if exist(nameseq) == 0
               file = nameseq;
            else
               error('Sequence cannot be found!');
            end
         end
         this.k = 1;
         this.sz = numel(dir([this.imgfolder '/*.png']));
         if this.sz > 0
            this.mask = ['%0' int2str(this.digits) 'd.png'];
         else
            this.sz = numel(dir([this.imgfolder '/*.jpg']));
            if this.sz > 0
               this.mask = ['%0' int2str(this.digits) 'd.jpg'];
            else
               this.sz = numel(dir([this.imgfolder '/*.tiff']));
               if this.sz > 0
                  this.mask = ['%0' int2str(this.digits) 'd.tiff'];
               else
                  error('Unknown file ending for sequence frames');
               end
            end
         end
         this.mask = [this.imgfolder '/' this.mask];
         frm = this.get_frame(1);
         [this.h,this.w,~] = size(frm);
      end

      function frame = get_next(this)
         frame = [];
         if this.has_next()
            frame = this.get_frame(this.k);
            this.k = this.k + 1;
         end
      end

      function has = has_next(this)
         has = this.k <= this.size();
      end

      function [framed, frame] = get_frame(this, n)
         frame = imread(sprintf(this.mask, n));
         framed = im2double(frame);
         if ~isempty(this.apply_function)
            frame = this.apply_function(frame);
            framed = this.apply_function(framed);
         end
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

