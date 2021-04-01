module SavGol
using DSP, LinearAlgebra

export SG,apply_filter

function vandermonde(halfwindow::Int, polyDeg::Int, T::Type=Float64)
    @assert halfwindow>=0
    @assert polyDeg>=0
    x = T[i for i in -halfwindow:halfwindow]
    n = polyDeg+1
    m = length(x)
    V = Array{T}(undef, m,n)
    for i in 1:m
        V[i,1]=T(1)
    end
    for j in 2:n
        for i in 1:m
            V[i,j] = x[i] *V[i,j-1]
        end
    end
    return V
end

function SG(halfWindow::Int, polyDeg::Int, T::Type=Float64)
    @assert 2*halfWindow > polyDeg

    V = vandermonde(halfWindow, polyDeg,T)
    Q,R = qr(V)
    D = [factorial(i) for i in 0:polyDeg] |> Diagonal
    SG = D*R \ Matrix(I, size(V')...) * Q'
    return SG'
end

function apply_filter(filter::Vector, signal::Vector)
    @assert isodd(length(filter))

    halfWindow = round(Int, (length(filter)-1)/2)
    padded_signal = vcat(signal[1]*ones(halfWindow),signal,signal[end] * ones(halfWindow))
    filter_cross_signal = conv(filter[end:-1:1], padded_signal)
    return filter_cross_signal[2*halfWindow+1:end-2*halfWindow]
end

end # module
