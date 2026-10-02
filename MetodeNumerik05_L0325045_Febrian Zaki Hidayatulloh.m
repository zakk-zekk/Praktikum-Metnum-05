A = [2 1 -1;
     4 3 1;
     -2 1 2];

b = [3; 9; 4];

n = length(b);

fprintf('========================================\n');
fprintf('PENYELESAIAN SPL DENGAN 3 METODE\n');
fprintf('========================================\n\n');

% =========================================================
% 1. ELIMINASI GAUSS
% =========================================================

Aug = [A b];

for k = 1:n-1
    for i = k+1:n
        m = Aug(i,k) / Aug(k,k);
        Aug(i,:) = Aug(i,:) - m * Aug(k,:);
    end
end

x_gauss = zeros(n,1);

for i = n:-1:1
    x_gauss(i) = (Aug(i,n+1) - Aug(i,i+1:n) * x_gauss(i+1:n)) / Aug(i,i);
end

fprintf('1. ELIMINASI GAUSS\n');
fprintf('Matriks Segitiga Atas:\n');
disp(Aug);

fprintf('Hasil:\n');
fprintf('x1 = %.4f\n', x_gauss(1));
fprintf('x2 = %.4f\n', x_gauss(2));
fprintf('x3 = %.4f\n\n', x_gauss(3));


% =========================================================
% 2. GAUSS-JORDAN
% =========================================================

Aug = [A b];

for k = 1:n
    Aug(k,:) = Aug(k,:) / Aug(k,k);

    for i = 1:n
        if i ~= k
            m = Aug(i,k);
            Aug(i,:) = Aug(i,:) - m * Aug(k,:);
        end
    end
end

x_gj = Aug(:,n+1);

fprintf('2. GAUSS-JORDAN\n');
fprintf('Matriks Identitas:\n');
disp(Aug);

fprintf('Hasil:\n');
fprintf('x1 = %.4f\n', x_gj(1));
fprintf('x2 = %.4f\n', x_gj(2));
fprintf('x3 = %.4f\n\n', x_gj(3));


% =========================================================
% 3. DEKOMPOSISI LU
% =========================================================

L = eye(n);
U = A;

for k = 1:n-1
    for i = k+1:n
        m = U(i,k) / U(k,k);
        L(i,k) = m;
        U(i,:) = U(i,:) - m * U(k,:);
    end
end

% Forward Substitution: Ly = b
y = zeros(n,1);

for i = 1:n
    y(i) = b(i) - L(i,1:i-1) * y(1:i-1);
end

% Backward Substitution: Ux = y
x_lu = zeros(n,1);

for i = n:-1:1
    x_lu(i) = (y(i) - U(i,i+1:n) * x_lu(i+1:n)) / U(i,i);
end

fprintf('3. DEKOMPOSISI LU\n');

fprintf('Matriks L:\n');
disp(L);

fprintf('Matriks U:\n');
disp(U);

fprintf('Nilai y:\n');
disp(y);

fprintf('Hasil:\n');
fprintf('x1 = %.4f\n', x_lu(1));
fprintf('x2 = %.4f\n', x_lu(2));
fprintf('x3 = %.4f\n\n', x_lu(3));


% =========================================================
% PERBANDINGAN HASIL
% =========================================================

fprintf('========================================\n');
fprintf('PERBANDINGAN HASIL\n');
fprintf('========================================\n');

fprintf('             x1        x2        x3\n');
fprintf('Gauss      %.4f    %.4f    %.4f\n', ...
        x_gauss(1), x_gauss(2), x_gauss(3));

fprintf('Gauss-Jordan %.4f    %.4f    %.4f\n', ...
        x_gj(1), x_gj(2), x_gj(3));

fprintf('LU         %.4f    %.4f    %.4f\n', ...
        x_lu(1), x_lu(2), x_lu(3));
