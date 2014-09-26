from http://users.rowan.edu/~polikar/WAVELETS/WTpart4.html:

the time localization will have a resolution that depends on which level they appear. If the main information of the signal lies in the high frequencies, as happens most often, the time localization of these frequencies will be more precise, since they are characterized by more number of samples. If the main information lies only at very low frequencies, the time localization will not be very precise, since few samples are used to express signal at these frequencies. This procedure in effect offers a good time resolution at high frequencies, and good frequency resolution at low frequencies. Most practical signals encountered are of this type.

https://www.mathworks.com/matlabcentral/newsreader/view_thread/267652

Process:
-generate a noisy signal to analyze
-test in matlab or python different WTs or DFTs as per link^
-meanwhile/ after consider programming and memory duress
-can WT run continuously while discarding old data?
-alternative run DFT on overlapping windows