function [y,dy] = trun_da(x,s)
% these are only here for trial_generator. find_intersect_truncated_cats uses an inline version for speed. this is redundant and annoying

%  maybe don't pre-define the functions, since they're only used once.

% these guys are pre-calculated and passed in a struct. just listing here for reference.
% s.t=sqrt(log(sig2/sig1)/(.5*(1/sig1^2 - 1/sig2^2)));
% s.sigsq = sig^2;
% s.sigsqrt2=s.sigsq*sqrt(2);
% s.sigsq1 = s.sigsq + sig1^2;
% s.sigsq2 = s.sigsq + sig2^2;
% s.sig1sqrt2 = sig1*sqrt(2);
% s.sig2sqrt2 = sig2*sqrt(2);
% s.k_1 = sig.*sig1./sqrt(s.sigsq1);
% s.k_2 = sig.*sig2./sqrt(s.sigsq2);
% s.k1sq = s.k_1.^2;
% s.k2sq = s.k_2.^2;
% s.tk_1sqrt2=s.t/(s.k_1*sqrt(2));
% s.tk_2sqrt2=s.t/(s.k_2*sqrt(2));
% s.k1s2 = s.k_1./s.sigsqrt2;
% s.k2s2 = s.k_2./s.sigsqrt2;

% d(x)
d = @(x) log(sqrt(s.sigsq2/s.sigsq1)) ...
    + log(1-erf(s.t/s.sig2sqrt2)) ...
    - log(erf(s.t./s.sig1sqrt2)) ...
    + (.5*x.^2)*(1/s.sigsq2-1/s.sigsq1) ...
    + log(erf(s.tk_1sqrt2 + x*s.k1s2) + erf(s.tk_1sqrt2 - x*s.k1s2)) ...
    - log(2 - erf(s.tk_2sqrt2 + x*s.k2s2) - erf(s.tk_2sqrt2 - x*s.k2s2));

% approximation for large pos/neg x
da = @(x) log(sqrt(s.sigsq2/s.sigsq1)) ...
    + log(1-erf(s.t/s.sig2sqrt2)) ...
    - log(erf(s.t./s.sig1sqrt2)) ...
    + (.5*x.^2)*(1/s.sigsq2-1/s.sigsq1) ...
    - s.t^2/(2*s.k1sq) - x.^2 * s.k1sq / (2*s.sig^4) ...
    + sign(x)*s.t.*x/s.sigsq ... % sign change here
    -log(-s.tk_1sqrt2 +sign(x).*x*s.k1s2)... % and here.
    -log(2*sqrt(pi));% slightly faster than -.5*log(pi)-log(2)

% approximation for x->0. this is a bit off.
da0 = @(x) log(sqrt(s.sigsq2/s.sigsq1)) ...
    + log(1-erf(s.t/s.sig2sqrt2)) ...
    - log(erf(s.t./s.sig1sqrt2)) ...
    + (.5*x.^2)*(1/s.sigsq2-1/s.sigsq1) ...
    + log(erf(s.tk_1sqrt2 + x*s.k1s2) + erf(s.tk_1sqrt2 - x*s.k1s2)) ...
    + .5 * (log(pi)-log(2)) ...
    - log(s.k_2) ...
    - 2*log(s.sig) ...
    + s.t^2 / (2*s.k2sq) ...
    + x.^2 * s.k2sq ./ (2*s.sig^4) ...
    - log(exp(-s.t*x/s.sigsq)./(s.t*s.sigsq - x.*s.k2sq)+exp(s.t.*x./s.sigsq)./(s.t*s.sigsq + x.*s.k2sq));

% with erfaterms
% da0 = @(x) log(sqrt(s.sigsq2/s.sigsq1)) ...
%     + log(1-erf(s.t/s.sig2sqrt2)) ...
%     - log(erf(s.t./s.sig1sqrt2)) ...
%     + (.5*x.^2)*(1/s.sigsq2-1/s.sigsq1) ...
%     + log(erf(s.tk_1sqrt2 + x*s.k1s2) + erf(s.tk_1sqrt2 - x*s.k1s2)) ...
%     + .5 * (log(pi)-log(2)) ...
%     - log(s.k_2) ...
%     - 2*log(s.sig) ... %2x removed but not for any good reason
%     + s.t^2 / (2*s.k2sq) ...
%     + x.^2 * s.k2sq ./ (2*s.sig^4) ...
%     - log(erfaterms(s.t.*s.sigsq + x.*s.k2sq,999).*exp(-s.t*x/s.sigsq)./(s.t*s.sigsq - x.*s.k2sq)+erfaterms(s.t.*s.sigsq - x.*s.k2sq,999).*exp(s.t.*x./s.sigsq)./(s.t*s.sigsq + x.*s.k2sq));


y = zeros(size(x));

erf_plus_erf=erf(s.tk_1sqrt2 + x*s.k1s2) + erf(s.tk_1sqrt2 - x*s.k1s2);
two_minus_erf_minus_erf = 2 - erf(s.tk_2sqrt2 + x*s.k2s2) - erf(s.tk_2sqrt2 - x*s.k2s2);

tol=1e-10;
idx_in=find(erf_plus_erf > tol & two_minus_erf_minus_erf > tol);
x_in = x(idx_in);
y(idx_in) = d(x_in);

idx_approx_large=find(erf_plus_erf <= tol);
x_out_large = x(idx_approx_large);
y(idx_approx_large) = da(x_out_large);

idx_approx_small=find(two_minus_erf_minus_erf <= tol);
x_out_small = x(idx_approx_small);
y(idx_approx_small) = da0(x_out_small);



if nargout==2
    % d/dx d(x)
    dd = @(x) (2/(s.sigsq*sqrt(2*pi))) * ( ...
        s.k_1*(exp(-(s.tk_1sqrt2 + x*s.k1s2).^2) - exp(-(s.tk_1sqrt2 - x*s.k1s2).^2)) ./ (erf(s.tk_1sqrt2 + x*s.k_1/s.sigsqrt2) + erf(s.tk_1sqrt2 - x*s.k_1/s.sigsqrt2)) ...
        + s.k_2*(exp(-(s.tk_2sqrt2 + x*s.k2s2).^2) - exp(-(s.tk_2sqrt2 - x*s.k2s2).^2)) ./ (2 - erf(s.tk_2sqrt2 + x*s.k2s2) - erf(s.tk_2sqrt2 - x*s.k2s2))) ...
        + x*(1/s.sigsq2-1/s.sigsq1);
    
    % approximation for large pos/neg x
    dda= @(x) x*(1/s.sigsq2-1/s.sigsq1) ...
        - x*s.k1sq/s.sig^4 ...
        + sign(x) * s.t/s.sigsq ...
        - s.k1sq ./ (-s.t*s.sigsq + x*s.k1sq);
    
    % approximation for x->0. this is a bit off.
    dda0 = @(x) (2/(s.sigsq*sqrt(2*pi))) * ( ...
        s.k_1*(exp(-(s.tk_1sqrt2 + x*s.k1s2).^2) - exp(-(s.tk_1sqrt2 - x*s.k1s2).^2)) ./ (erf(s.tk_1sqrt2 + x*s.k_1/s.sigsqrt2) + erf(s.tk_1sqrt2 - x*s.k_1/s.sigsqrt2)) ...
        + sqrt(pi/2)*(exp(-(s.tk_2sqrt2 + x*s.k2s2).^2) - exp(-(s.tk_2sqrt2 - x*s.k2s2).^2)) ...
        ./ (s.sigsq.*exp(-(s.t^2./(2*s.k2sq)) - x.^2 .* s.k2sq ./(2.*s.sig^4)) .* (exp(-s.t*x./s.sigsq)./(s.t.*s.sigsq + x.*s.k2sq) + exp(s.t*x./s.sigsq)./(s.t.*s.sigsq - x.*s.k2sq)))) ...
        + x*(1/s.sigsq2-1/s.sigsq1);
    
    % With erfa terms
    % dda0 = @(x) (2/(s.sigsq*sqrt(2*pi))) * ( ...
    %       s.k_1*(exp(-(s.tk_1sqrt2 + x*s.k1s2).^2) - exp(-(s.tk_1sqrt2 - x*s.k1s2).^2)) ./ (erf(s.tk_1sqrt2 + x*s.k_1/s.sigsqrt2) + erf(s.tk_1sqrt2 - x*s.k_1/s.sigsqrt2)) ...
    %     + sqrt(pi/2)*(exp(-(s.tk_2sqrt2 + x*s.k2s2).^2) - exp(-(s.tk_2sqrt2 - x*s.k2s2).^2)) ...
    %     ./ (s.sigsq.*exp(-(s.t^2./(2*s.k2sq)) - x.^2 .* s.k2sq ./(2.*s.sig^4)) .* (erfaterms(s.t.*s.sigsq + x.*s.k2sq,999).*exp(-s.t*x./s.sigsq)./(s.t.*s.sigsq + x.*s.k2sq) + erfaterms(s.t.*s.sigsq - x.*s.k2sq,999).*exp(s.t*x./s.sigsq)./(s.t.*s.sigsq - x.*s.k2sq)))) ...
    %     + x*(1/s.sigsq2-1/s.sigsq1);
    dy= zeros(size(x));
    dy(idx_in)= dd(x_in);
    dy(idx_approx_large)= dda(x_out_large);
    dy(idx_approx_small)= dda0(x_out_small);
end