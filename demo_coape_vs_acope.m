%%
% Why do we need to rectify:
% - there is an ambiguity between the dipole direction and the sign of the
% reconstructed time series
% - not trivial to resolve this
% - hence we need to find a way to do tests/comparisons that are
% insensitive to this ambiguity so that the statistic/copes are smooth in 
% space AND comparable over subjects.
%
% Options are:
% acope : uses abs(c'*pe), and needs to subtract the baseline average
% abs(cope) even for differential contrasts
% coape : uses c'*abs(pe), and needs to subtract the baseline average
% abs(cope) only for average (main effect) contrasts
% acope always does OK
% coape fails completely when one +ve and one -ve
% coape better than acope if both in same direction (e.g. both +ve)

% only one active
x=2+randn(100000,1);y=randn(100000,1);z=4+randn(100000,1);

figure; subplot(1,3,1);
[h,t]=hist(abs(x)-abs(y),50);[h2]=hist(abs(z)-abs(y),t);plot(t,h,'LineWidth',2);ho; plot(t,h2,'r','LineWidth',2); 
tpr_at_fpr05=sum(abs(z)-abs(y)>percentile(abs(x)-abs(y),95))/length(abs(z)-abs(y));
title(['coape: abs(x)-abs(y), ', num2str(tpr_at_fpr05)]);

subplot(1,3,2);
[h,t]=hist(abs(x-y),50);[h2]=hist(abs(z-y),t);plot(t,h,'LineWidth',2);ho; plot(t,h2,'r','LineWidth',2); 
tpr_at_fpr05=sum(abs(z-y)>percentile(abs(x-y),95))/length(abs(z-y));
title(['acope: abs(x-y), ', num2str(tpr_at_fpr05)]);

subplot(1,3,3);
[h,t]=hist(x-y,50);[h2]=hist(z-y,t);plot(t,h,'LineWidth',2);ho; plot(t,h2,'r','LineWidth',2); 
tpr_at_fpr05=sum((z-y)>percentile((x-y),95))/length((z-y));
title(['cope: x-y, ', num2str(tpr_at_fpr05)]);

% both active
x=2+randn(100000,1);y=2*randn(100000,1);z=4+randn(100000,1);

figure; subplot(1,3,1);
[h,t]=hist(abs(x)-abs(y),50);[h2]=hist(abs(z)-abs(y),t);plot(t,h,'LineWidth',2);ho; plot(t,h2,'r','LineWidth',2); 
tpr_at_fpr05=sum(abs(z)-abs(y)>percentile(abs(x)-abs(y),95))/length(abs(z)-abs(y));
title(['coape: abs(x)-abs(y), ', num2str(tpr_at_fpr05)]);

subplot(1,3,2);
[h,t]=hist(abs(x-y),50);[h2]=hist(abs(z-y),t);plot(t,h,'LineWidth',2);ho; plot(t,h2,'r','LineWidth',2); 
tpr_at_fpr05=sum(abs(z-y)>percentile(abs(x-y),95))/length(abs(z-y));
title(['acope: abs(x-y), ', num2str(tpr_at_fpr05)]);

subplot(1,3,3);
[h,t]=hist(x-y,50);[h2]=hist(z-y,t);plot(t,h,'LineWidth',2);ho; plot(t,h2,'r','LineWidth',2); 
tpr_at_fpr05=sum((z-y)>percentile((x-y),95))/length((z-y));
title(['cope: x-y, ', num2str(tpr_at_fpr05)]);

% one positive, one negative
x=-2+randn(100000,1);y=-2+randn(100000,1);z=2+randn(100000,1);

figure; subplot(1,3,1);
[h,t]=hist(abs(x)-abs(y),50);[h2]=hist(abs(z)-abs(y),t);plot(t,h,'LineWidth',2);ho; plot(t,h2,'r','LineWidth',2); 
tpr_at_fpr05=sum(abs(z)-abs(y)>percentile(abs(x)-abs(y),95))/length(abs(z)-abs(y));
title(['coape: abs(x)-abs(y), ', num2str(tpr_at_fpr05)]);

subplot(1,3,2);
[h,t]=hist(abs(x-y),50);[h2]=hist(abs(z-y),t);plot(t,h,'LineWidth',2);ho; plot(t,h2,'r','LineWidth',2); 
tpr_at_fpr05=sum(abs(z-y)>percentile(abs(x-y),95))/length(abs(z-y));
title(['acope: abs(x-y), ', num2str(tpr_at_fpr05)]);

subplot(1,3,3);
[h,t]=hist(x-y,50);[h2]=hist(z-y,t);plot(t,h,'LineWidth',2);ho; plot(t,h2,'r','LineWidth',2); 
tpr_at_fpr05=sum((z-y)>percentile((x-y),95))/length((z-y));
title(['cope: x-y, ', num2str(tpr_at_fpr05)]);


if(0)

%%%% using baseline as null
%%


% both active positive
w=randn(100000,1);x=randn(100000,1);y=2+randn(100000,1);z=4+randn(100000,1);

[h,t]=hist(abs(x)-abs(w),50);[h2]=hist(abs(z)-abs(y),t);figure;plot(t,h);ho; plot(t,h2,'r'); 
tpr_at_fpr05=sum(abs(z)-abs(y)>percentile(abs(x)-abs(w),95))/length(abs(z)-abs(y));
title(['coape: abs(x)-abs(y), ', num2str(tpr_at_fpr05)]);
ca

[h,t]=hist(abs(x-w),50);[h2]=hist(abs(z-y),t);figure;plot(t,h);ho; plot(t,h2,'r'); 
tpr_at_fpr05=sum(abs(z-y)>percentile(abs(x-w),95))/length(abs(z-y));
title(['acope: abs(x-y), ', num2str(tpr_at_fpr05)]);



% only one active
w=randn(100000,1);x=randn(100000,1);y=randn(100000,1);z=2+randn(100000,1);

[h,t]=hist(abs(x)-abs(w),50);[h2]=hist(abs(z)-abs(y),t);figure;plot(t,h);ho; plot(t,h2,'r'); 
tpr_at_fpr05=sum(abs(z)-abs(y)>percentile(abs(x)-abs(w),95))/length(abs(z)-abs(y));
title(['coape: abs(x)-abs(y), ', num2str(tpr_at_fpr05)]);

[h,t]=hist(abs(x-w),50);[h2]=hist(abs(z-y),t);figure;plot(t,h);ho; plot(t,h2,'r'); 
tpr_at_fpr05=sum(abs(z-y)>percentile(abs(x-w),95))/length(abs(z-y));
title(['acope: abs(x-y), ', num2str(tpr_at_fpr05)]);


% one positive, one negative
w=randn(100000,1);x=randn(100000,1);y=-2+randn(100000,1);z=2+randn(100000,1);

[h,t]=hist(abs(x)-abs(w),50);[h2]=hist(abs(z)-abs(y),t);figure;plot(t,h);ho; plot(t,h2,'r'); 
tpr_at_fpr05=sum(abs(z)-abs(y)>percentile(abs(x)-abs(w),95))/length(abs(z)-abs(y));
title(['coape: abs(x)-abs(y), ', num2str(tpr_at_fpr05)]);

[h,t]=hist(abs(x-w),50);[h2]=hist(abs(z-y),t);figure;plot(t,h);ho; plot(t,h2,'r'); 
tpr_at_fpr05=sum(abs(z-y)>percentile(abs(x-w),95))/length(abs(z-y));
title(['acope: abs(x-y), ', num2str(tpr_at_fpr05)]);


% null
w=randn(100000,1);x=randn(100000,1);y=4+randn(100000,1);z=4+randn(100000,1);

[h,t]=hist(abs(x)-abs(w),50);[h2]=hist(abs(z)-abs(y),t);figure;plot(t,h);ho; plot(t,h2,'r'); 
tpr_at_fpr05=sum(abs(z)-abs(y)>percentile(abs(x)-abs(w),95))/length(abs(z)-abs(y));
title(['coape: abs(x)-abs(y), ', num2str(tpr_at_fpr05)]);

[h,t]=hist(abs(x-w),50);[h2]=hist(abs(z-y),t);figure;plot(t,h);ho; plot(t,h2,'r'); 
tpr_at_fpr05=sum(abs(z-y)>percentile(abs(x-w),95))/length(abs(z-y));
title(['acope: abs(x-y), ', num2str(tpr_at_fpr05)]);

end;
