%cycles through the components of DSTREAM... and displays them one after
%another. Input: DATA . Optional input: PHASE


function bss_disp(DATA,varargin)

if ~isempty(varargin)
    PHASE = varargin{1};
end

%fprintf(['plotting data in: ' str(DATA)])

[s1,s2] = size(DATA);

if s1 < s2 
    DATA = transpose(DATA);
end

[s1,s2] = size(DATA);

i=1;
while i ~= (s2+1)
  GOON = input('Next component or specific? y/n/NUM [y]: ','s');
  if isempty(GOON)
      GOON = 'y';
  elseif GOON == 'n'
       break
  end
  if GOON =='y'
      figure(20) 
      clf()
      if ~isempty(varargin)
          plot(PHASE,DATA(:,i),'x')
      else
          plot(DATA(:,i),'x')
      end
      title(['component: ' num2str(i)])
      fprintf(['component: ' num2str(i) '\n']) 
  else
      figure(20) 
      clf()
      if ~isempty(varargin)
          plot(PHASE,DATA(:,str2num(GOON)),'x')
      else
          plot(DATA(:,str2num(GOON)),'x')
      end
      title(['component: ' GOON])
      fprintf(['component: ' GOON '\n']) 
      i = str2num(GOON);
  end
  i = i+1;
end



