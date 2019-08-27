function h = flipbits(bitmask, g)
    bitmask = int64(bitmask);
    bits(5) = mod(idivide(bitmask, int64(16)), 2);
    bits(4) = mod(idivide(bitmask, int64(8)), 2);
    bits(3) = mod(idivide(bitmask, int64(4)), 2);
    bits(2) = mod(idivide(bitmask, int64(2)), 2);
    bits(1) = mod(bitmask, 2); 
    syms b_1 b_2 b_3 b_4 b_a k x y z w
    vars = [b_1, b_2, b_3, b_4];
    coeffs = [y, z, w, x, k];
    foriginal = k*b_1*b_2*b_3*b_4 + x*b_1*b_2*b_3 + y*b_2*b_3*b_4 + z*b_3*b_4*b_1 + w*b_4*b_1*b_2;
    f = foriginal;
    if bits(5) == 1
        g = subs(g, {b_a}, {1-b_a});
    end
    for i = 1:4 
        if bits(i) == 1
            g = subs(g, {vars(i)}, {1-vars(i)});
            f = subs(f, {vars(i)}, {1-vars(i)});

            for j = 1:5
                if j ~= i 
                    g = subs(g, {coeffs(j)}, {-coeffs(j)});
                    f = subs(f, {coeffs(j)}, {-coeffs(j)});
                end
            end
            g = subs(g, {coeffs(i)}, {k + coeffs(i)});
            f = subs(f, {coeffs(i)}, {k + coeffs(i)});
        end
    end
    g = expand(g- f + foriginal);
    g = collect(g, k);
    g = collect(g, x);
    g = collect(g, y);
    g = collect(g, z);
    g = collect(g, w);
    h = collect(g, b_a);
end
