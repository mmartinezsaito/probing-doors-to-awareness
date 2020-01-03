function plt_spectra(T, D)

  ns = length(T);
  TF = linspace(0, 1, ns/2) * nyl;
  DF = fft(D);
  DF = DF(1:ns/2);
  psdD = abs(DF).^2 / (fs*ns);
  psdD(2:end-1) = 2*psdD(2:end-1);
  phi = atan(imag(DF)./real(DF));
  
  % periodogram PSD and phase
  fh = figure;  
  ah(1) = subplot(2, 3, 1); plot(TF, psdD), axis([0 2 0 10])  % power spectrum
  ah(2) = subplot(2, 3, 4); semilogy(TF, psdD) 
  ah(3) = subplot(2, 3, 2); plot(TF, 10*log10(psdD)), title('Periodogram using FFT'), grid on
                            xlabel('Frequency (Hz)'), ylabel('Power/Frequency (dB/Hz)')
  ah(3) = subplot(2, 3, 5); periodogram(DS, rectwin(ns), ns, fs, 'psd')
  ah(4) = subplot(2, 3, 3); plot(TF, phi(1:ns/2));
  ah(3) = subplot(2, 3, 6); semilogy(TF, phi(1:ns/2)); 
  
  % PSD estimate using Welch's averaged, modified periodogram method
  fh = figure;  
  Pxx = pwelch(D, [], [], [], fs, 'centered');
  plot(Pxx)

  % spectrogram
  fh = figure; 
  ah(1) = subplot(2,2,1); 
  spectrogram(D0, 256, [], [], fs, 'yaxis'), colorbar
  ah(2) = subplot(2,2,2); 
  spectrogram(D1, 256, [], [], fs, 'yaxis'), colorbar
  ah(3) = subplot(2,2,3); 
  spectrogram(D2, 256, [], [], fs, 'yaxis'), colorbar
  ah(4) = subplot(2,2,4);
  spectrogram(DS, 256, [], [], fs, 'yaxis'), colorbar
