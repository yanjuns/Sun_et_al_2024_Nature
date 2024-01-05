%t-time values over which we wish to calculate the derivative. The time
%values should be equidistant, or as close as possible.

%P-Position Data

%order should be a relatively low degree polynomial to prevent overfitting
%framelen needs to be odd, 
%I have used order=5, framelen=21 in the test script

%This method seems to have issues near the endpoints, at least for the
%function I tried.
function y=sgolayDerivative(t,P,order,framelen)
%Use the average distance between t values, if they  aren't exactly
%equidistant. If they are, this will yield the distance between each point
dt=mean(t(2:end)-t(1:end-1))
%Apply the savgol filter
[b,g]=sgolay(order,framelen);
D=zeros(length(t),2);
% size(P)
% size(g(:,0+1))
for p=0:1
    D(:,p+1)=conv(P,factorial(p)/(-dt)^p*g(:,p+1),'same');
end

%Comment these out if you don't want to see the plots
%figure()
%plot(t,P,t,D(:,1))
%legend('Position','Smoothed Position')
%figure()
%plot(t,D(:,2))
%legend('Approximate Derivative')

y=D(:,2);