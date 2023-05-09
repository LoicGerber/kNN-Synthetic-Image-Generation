function spae = spae_metric(im1, im2)
    % Compute the phase-only matched filter (POMF) images of the two input images
    pomf1 = real(ifft2(abs(fft2(im1)).*exp(1i*angle(fft2(im2)))));
    pomf2 = real(ifft2(abs(fft2(im2)).*exp(1i*angle(fft2(im1)))));

    % Compute the absolute error between the POMF images
    spae = sum(abs(pomf1(:)-pomf2(:))) / numel(pomf1);
end
