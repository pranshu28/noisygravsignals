using HDF5
using PyPlot
using DSP
using StatsBase
#using PyCall
#@pyimport numpy as np

#H1_16_32 = h5open("/Data/H-H1_LOSC_16_V1-1126259446-32.hdf5", "r");
H1_4_32 = h5open("/Data/H-H1_LOSC_4_V1-1126259446-32.hdf5", "r");
#L1_16_32 = h5open("/Data/L-L1_LOSC_16_V1-1126259446-32.hdf5", "r");
L1_4_32 = h5open("/Data/L-L1_LOSC_4_V1-1126259446-32.hdf5", "r");
GW=readdlm("/Data/GW150914_4_NR_waveform.txt");

n=10000
#s1=read(H1_16_32["strain"])["Strain"]
s2=read(H1_4_32["strain"])["Strain"][1:n]
#s3=read(L1_16_32["strain"])["Strain"]
s4=read(L1_4_32["strain"])["Strain"][1:n]
m2=read(H1_4_32["meta"])["Duration"]
t2=read(H1_4_32["meta"])["GPSstart"]

#close(H1_16_32)
close(H1_4_32);
#close(L1_16_32)
close(L1_4_32);

#Plot the detector strains and the template for comparison
nt = length(s2)
interval = m2/nt
xlabel("Time")
ylabel("Strain")
title("Raw")
#x=collect(t2:interval:t2+m2-interval)
#plot(x,s2)
plot(s4)
#plot(GW[:,1],GW[:,2])


#Removing the noise

nft = round(Int,1/interval)
h_psd,freq1 = psd(s2,nft,nft)         #psd1
l_psd,freq2 = psd(s4,nft,nft)         #psd2
temp_psd,freq3 = psd(GW[:,2],nft,nft) #temp_psd
loglog(freq1,sqrt(h_psd))  #asd1
loglog(freq2,sqrt(l_psd))  #asd2
loglog(freq3,sqrt(temp_psd))  #temp_asd
#Whitening
function whiten(strain,i)
  ft = fft(strain)
  pd,freq = psd(strain,nft,nft)
  white_ft = ft / (sqrt(pd /i/2.))                  #####################
  white_strain = ifft(white_ft)
  return white_strain
end
hstrain_whiten = whiten(s2, interval)
lstrain_whiten = whiten(s4, interval)
temp_whiten = whiten(GW[:,2], interval)

responsetype = Bandpass(20/(nft/2.0), 300/(nft/2.0),fs=1000)
designmethod = Butterworth(4)
h_red = filt(digitalfilter(responsetype, designmethod), hstrain_whiten)
l_red = filt(digitalfilter(responsetype, designmethod), lstrain_whiten)
temp_red = filt(digitalfilter(responsetype, designmethod), temp_whiten)

#Plot Whitened data
xlabel("Time")
ylabel("Strain")
title("Whitened")
#plot(h_red)
plot(l_red)
#plot(temp_red)


#Finding the signal
