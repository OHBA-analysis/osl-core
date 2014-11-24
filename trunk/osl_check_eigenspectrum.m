function dim2use = osl_check_eigenspectrum(allsvd,pcadim,do_plots)

%% Check the eigenspectrum for discontinuities introduced by MaxFilter and ICA preprocessing

diff_svd=diff(log(allsvd));
diff_svd_n = abs(diff_svd/log(allsvd(1)));

kicks=find(diff_svd_n > mean(diff_svd_n)+5*std(diff_svd_n));

if ~isempty(kicks), kicks = kicks(kicks>40); end
for nk=1:numel(kicks)
    tmp=abs(diff_svd(kicks(nk)-5:kicks(nk)));
    [~,tmp2]=max(tmp);
    kicks(nk)=kicks(nk)-6+tmp2;
end  
if ~isempty(kicks) && kicks(1) < pcadim
    if do_plots; figure; plot(log(allsvd)); ho; plot(kicks, log(allsvd(kicks)), 'rx');ho; plot(pcadim,log(allsvd(pcadim)), 'go'); title('Eigenspectrum and rank of pseudo-inverse');legend('Eigenspectrum','Location of kicks','Original rank-estimate'); end
    fprintf('\n%s%d%s%d%s\n', 'The estimated dimensionality has been forced from ',pcadim,' down to ',kicks(1),'.');
    dim2use = kicks(1);
else
    if do_plots; figure; plot(log(allsvd)); ho; plot(pcadim,log(allsvd(pcadim)), 'go'); title('Eigenspectrum and rank of pseudo-inverse');legend('Eigenspectrum','Original rank-estimate'); end
    fprintf('\n%s%d%s\n', 'No kicks found in the eigenspectrum. Estimated dimensionality has been kept at ',pcadim,'.');
    dim2use = pcadim;
end
end