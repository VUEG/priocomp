\ Model multiobj
\ LP format - for model browsing. Use MPS format to capture full model detail.
Maximize multi-objectives
  Value: Priority=0 Weight=1 AbsTol=1 RelTol=0.01
   El[0] + 4 El[1] + El[2] + El[3] + El[4] + El[5] + El[6] + El[7] + El[8]
   + 4 El[9] + 4 El[13] + 4 El[16] + 4 El[19]
  Cost: Priority=0 Weight=-1 AbsTol=2 RelTol=0.01
   2 El[1] + El[5] + El[6] + El[7] + El[8] + 2 El[9] + 2 El[13] + El[15]
   + 3 El[16] + El[17] + El[18] + 2 El[19]
Subject To
 R0: 0.05 El[0] + 0.05 El[1] + 0.05 El[2] + 0.05 El[3] + 0.05 El[4]
   + 0.05 El[5] + 0.05 El[6] + 0.05 El[7] + 0.05 El[8] + 0.05 El[9]
   + 0.05 El[10] + 0.05 El[11] + 0.05 El[12] + 0.05 El[13] + 0.05 El[14]
   + 0.05 El[15] + 0.05 El[16] + 0.05 El[17] + 0.05 El[18] + 0.05 El[19]
   = 0.2
Bounds
Binaries
 El[0] El[1] El[2] El[3] El[4] El[5] El[6] El[7] El[8] El[9] El[10] El[11]
 El[12] El[13] El[14] El[15] El[16] El[17] El[18] El[19]
End
