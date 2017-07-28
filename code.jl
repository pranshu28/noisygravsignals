using HDF5
using PyPlot
using DSP
using Interpolations

#H1_16_32 = h5open("/dataset/H-H1_LOSC_16_V1-1126259446-32.hdf5", "r");
H1_4_32 = h5open("dataset/H-H1_LOSC_4_V1-1126259446-32.hdf5", "r");
#L1_16_32 = h5open("/dataset/L-L1_LOSC_16_V1-1126259446-32.hdf5", "r");
L1_4_32 = h5open("dataset/L-L1_LOSC_4_V1-1126259446-32.hdf5", "r");
GW=readdlm("dataset/GW150914_4_NR_waveform.txt");

#s1=read(H1_16_32["strain"])["Strain"]
s2=read(H1_4_32["strain"])["Strain"]
#s3=read(L1_16_32["strain"])["Strain"]
s4=read(L1_4_32["strain"])["Strain"]
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
#loglog(freq1,sqrt(h_psd))  #asd1
#loglog(freq2,sqrt(l_psd))  #asd2
#loglog(freq3,sqrt(temp_psd))  #temp_asd

h_itp = interpolate((freq1,), h_psd, Gridded(Linear()))
l_itp = interpolate((freq2,), l_psd, Gridded(Linear()))

#Whitening
function whiten(strain,psd_,i)
  strain_len = length(strain)
  freq = rfftfreq(strain_len,i)
  hf = rfft(strain)
  white_ft = Complex{Float64}[]
  asd = sqrt(psd_[freq] /i/2.)
  for i=1:length(hf)
    push!(white_ft,hf[i]/asd[i])
  end
  white_strain = irfft(white_ft,strain_len)
  return white_strain
end
hstrain_whiten = whiten(s2, h_itp, interval)
lstrain_whiten = whiten(s4, l_itp, interval)
temp_whiten = whiten(GW[:,2], h_itp, interval)

#Plot Whitened data
xlabel("Time")
ylabel("Strain")
title("Whitened")
responsetype = Bandpass(20/(nft/2.0), 300/(nft/2.0))
designmethod = Butterworth(4)
h_red = filt(digitalfilter(responsetype, designmethod), hstrain_whiten)
#plot(h_red)
l_red = filt(digitalfilter(responsetype, designmethod), lstrain_whiten)
#plot(l_red)
temp_red = filt(digitalfilter(responsetype, designmethod), temp_whiten)
#plot(temp_red)


#Finding the signal

#Plot the correlation between each detection strain and the reference template
hcorr = xcorr(h_red,temp_red)
lcorr = xcorr(l_red,temp_red)
plot(hcorr)
plot(lcorr)

#Plot the whitened H1 strain-template strain correlation between 15 and 17 seconds
startind = round(Int, 127000*(131072/262143))
endind = round(Int, 135000*(131072/262143))
hcorr = xcorr(h_red[startind:endind], temp_red)
plot(hcorr)
lcorr = xcorr(l_red[startind:endind], temp_red)
plot(lcorr)

#Zero in on the event
startind = 1400
endind = 6600
plot(hcorr[startind:endind])
plot(lcorr[startind:endind])

#Using the time of the event provided by LIGO, see if you found it
xlabel("Time")
ylabel("Strain")
title("Whitened")
plot(h_red[startind:endind])
plot(l_red[startind:endind])
plot(temp_red)
