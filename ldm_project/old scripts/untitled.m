t = 360*60;

times = repmat(300,t/300,1);

lights = repmat([1;0],(t*0.5)/300,1);

inten = repmat([5 5 75 75 225 225]',t/(300*6),1);

designFile = [times lights inten];

dlmwrite('lightSequence_5Min_VariedIntensity.txt',designFile)

