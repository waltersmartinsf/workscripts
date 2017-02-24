%bss wrapper for multicombi algorithm 


function [DSTREAM, PSTREAM] = bss_mcombi(DSTREAM,PSTREAM)

%transposing signal if needed
[PCAs1,PCAs2] = size(DSTREAM.DATA.pca);
if PCAs1 > PCAs2
   DATAtrans = transpose(DSTREAM.DATA.pca);
else
   DATAtrans = DSTREAM.DATA.pca;
end

%running mutlicombi
if PSTREAM.MCOMBI.arcoeff ~= 0;
  [W, nonseparablecomponents, Wefica, Wwasobi, ISRwa1, ISRef1, signals]= multicombi(DATAtrans,PSTREAM.MCOMBI.arcoeff);
else
  [W, nonseparablecomponents, Wefica, Wwasobi, ISRwa1, ISRef1, signals]= multicombi(DATAtrans);
end
  
%testing whether output is complex. if so warn and revert
sigreal = isreal(signals);
if sigreal ~= 1
    fprintf('\n')
    fprintf('Warning: output signals are complex. Retaining real parts only. \n')
    fprintf('...')
    DSTREAM.INFO.WARNINGS.mcombi2 = 'Warning: output signals are complex. Retaining real parts only.';
    signals = real(signals);
end

DSTREAM.DATA.MCOMBI.w = W;
DSTREAM.DATA.MCOMBI.nonsep = nonseparablecomponents;
DSTREAM.INFO.MCOMBI.Wefica = Wefica;
DSTREAM.INFO.MCOMBI.Wwasobi = Wwasobi;
DSTREAM.INFO.MCOMBI.ISRwa1 = ISRwa1;
DSTREAM.INFO.MCOMBI.ISRef1 = ISRef1;
DSTREAM.DATA.MCOMBI.signals = signals;

