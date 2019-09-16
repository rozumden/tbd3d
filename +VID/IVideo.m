% Abstract class for creating Video
classdef (Abstract) IVideo < handle
   methods (Abstract)
      has = has_next(this)
      frame = get_next(this)
      frame = get_frame(this, n)
      sz = length(this)
      h = height(this)
      w = width(this)
      close(this)
   end
   
   methods
      function [] = show_frame(this, n)
         imshow(this.get_frame(n));
      end
   end
end 