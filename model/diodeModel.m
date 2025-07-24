function JD = diodeModel(V, params, config)
    %DIODEMODEL Compute diode current using Lambert W analytic formula.
    %   JD = diodeModel(V, params, config) returns the current density for the
    %   voltage vector V. PARAMS = [J0, Rs, Rsh, k, J02]. The implementation
    %   solves the single diode equation with series and shunt resistances
    %   using the Lambert W function and then adds the tunnelling and non-ohmic
    %   contributions. The calculation is vectorised over V.

    J0  = params(1);
    Rs  = params(2);
    Rsh = params(3);
    k   = params(4);
    J02 = params(5);

    if Rs <= 0
        error('物理参数错误: Rs必须为正值 (当前值: %.6e)', Rs);
    end
  
    Vt = config.physics.kb * config.physics.T / config.physics.q;
    n  = config.physics.n;

    % Lambert W based solution for the base diode + shunt equation
    x = (Rs .* Rsh .* J0) ./ (n * Vt * (Rs + Rsh)) .* ...
        exp((Rsh .* (Rs .* J0 + V)) ./ (n * Vt * (Rs + Rsh)));

    if exist('lambertw', 'file')
        w = lambertw(x);
    else
        w = arrayfun(@approximateLambertW, x);
    end

    JD_base = (n * Vt ./ Rs) .* w - (Rsh .* J0 - V) ./ (Rs + Rsh);

    % Voltage drop across the diode after accounting for series resistance
    Vdrop = V - JD_base .* Rs;

    JD = JD_base + J02 .* (exp(config.physics.A2 .* Vdrop) - 1) + ...
        k .* (abs(Vdrop) .^ config.physics.m) .* sign(Vdrop);
end

function w = approximateLambertW(x)
    % Approximate Lambert W for environments without the symbolic toolbox
    if x < 0
        w = 0;
    elseif x == 0
        w = 0;
    elseif x < 1
        w = x * (1 - x + 1.5 * x^2 - 2.667 * x^3 + 5.208 * x^4);
    else
        if x < 3
            w = 0.5;  
        else
            w = log(x) - log(log(x));
        end
        for i = 1:10
            ew = exp(w);
            w_next = w - (w * ew - x) / (ew + w * ew);
            if abs(w_next - w) < 1e-10
                w = w_next;
                break;
            end
            w = w_next;
        end
    end
end
