using HDF5
using PyPlot
using DSP
#using StatsBase
using Interpolations
#using PyCall
#@pyimport numpy as np

#H1_16_32 = h5open("/Data/H-H1_LOSC_16_V1-1126259446-32.hdf5", "r");
H1_4_32 = h5open("Data/H-H1_LOSC_4_V1-1126259446-32.hdf5", "r");
#L1_16_32 = h5open("/Data/L-L1_LOSC_16_V1-1126259446-32.hdf5", "r");
L1_4_32 = h5open("Data/L-L1_LOSC_4_V1-1126259446-32.hdf5", "r");
GW=readdlm("Data/GW150914_4_NR_waveform.txt");

n=20000
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
#xlabel("Time")
#ylabel("Strain")
#title("Raw")
#x=collect(t2:interval:t2+m2-interval)
#plot(s2)
#plot(s4)
#plot(GW[:,1],GW[:,2])


#Removing the noise

nft = round(Int,1/interval)
h_psd,freq1 = psd(s2,nft,nft)         #psd1
l_psd,freq2 = psd(s4,nft,nft)         #psd2
temp_psd,freq3 = psd(GW[:,2],nft,nft) #temp_psd
loglog(freq1,sqrt(h_psd))  #asd1
loglog(freq2,sqrt(l_psd))  #asd2
loglog(freq3,sqrt(temp_psd))  #temp_asd

h_itp = interpolate((freq1,), h_psd, Gridded(Linear()))
l_itp = interpolate((freq2,), l_psd, Gridded(Linear()))

#Whitening
function whiten(strain,psd_,i)
  strain_len = length(strain)
  freq = rfftfreq(strain_len,i)
  hf = rfft(strain)
  white_ft = hf/(sqrt(psd_[freq] /i/2.))
  white_strain = ifft(white_ft)
  print(psd_[freq])
  return white_strain
end
hstrain_whiten = whiten(s2, h_itp, interval)
lstrain_whiten = whiten(s4, l_itp, interval)
temp_whiten = whiten(GW[:,2], temp_psd, interval)

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
#plot(l_red)
plot(temp_red)


#Finding the signal
hcorr = xcorr(h_red,l_red)
lcorr = xcorr(l_red,temp_red)
plot(hcorr)
plot(lcorr)

startind = #np.where(htime==(min(htime)+15))[0][0]
endind = #np.where(htime==(min(htime)+17))[0][0]
hcorr = np.xcorr(h_red[startind:endind], l_red)
plot(hcorr)

startind = #np.where(htime==(min(htime)+16.25))[0][0]
endind = #np.where(htime==(min(htime)+16.5))[0][0]
plot(htime[startind:endind], hcorr[startind:endind])
plot(htime[startind:endind], lcorr[startind:endind])
