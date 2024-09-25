function get_index_p(x,J,L)
    floor(Int64,(x>=0) ? mod(x / L[1] * J[1] + 1 ,J[1]) : -mod(-x / L[1] * J[1] ,J[1])+J+1 )
    end
    
function get_index_p(x::Vector,J,L)
        j = Vector{Int64}(undef,length(J))
        for i in 1:length(J)
            j[i] = get_index_p(x[i],J[i],L[i])
        end
        return j[:]
end

function n(ϕ,ϕ_t,Box_x,J)
        V = volume(Box_x)
        dx = differentials(Box_x,J)
        n = real(im*sum(ϕ.*conj.(ϕ_t).-conj.(ϕ).*ϕ_t))/2#*prod(dx)/2)
        return n, n/prod(J)
end



function f(k,f_pars)
    (k0,δ) = f_pars
    if abs(k'*k - k0^2) < δ^2
        return (1.0 + 0.0im)/(2π*(k0*2*δ)) # normalized so that the integral is 1
    else
        return 0.0 + 0.0im
    end
end

function h(k,h_pars)
    k0 = h_pars
    if abs(k'*k)< k0^2
        return (1.0 + 0.0im)/2*π/k0^2
    else
        return 0.0
    end
end

function get_frequencies(Box,J)
    κ = 2π./(Box_x[2:2:end] - Box_x[1:2:end-1])
    kx = fftfreq(J[1]).*κ[1]*J[1]
    ky = fftfreq(J[2]).*κ[2]*J[2]
    kz = fftfreq(J[3]).*κ[3]*J[3]
    return kx, ky, kz
end

ω(k,m2) = sqrt((k'*k+m2))

function get_fourier_data(f,Box_x,J,m2,f_pars)
    g = Array{ComplexF64}(undef,J...)
    g_p = copy(g)
    g_t = copy(g)
    κ = 2π./(Box_x[2:2:end] - Box_x[1:2:end-1])

    kx, ky, kz = get_frequencies(Box_x,J)
    
    kx = fftfreq(J[1]).*κ[1]*J[1]
    ky = fftfreq(J[2]).*κ[2]*J[2]
    kz = fftfreq(J[3]).*κ[3]*J[3]
    
    kx_p = sort(kx)
    ky_p = sort(ky)
    kz_p = sort(kz)

    for i in 1:J[1]
        for j in 1:J[2]
            for l in 1:J[3]
            ω_val = ω([kx[i];ky[j];kz[l]],m2)
            g[i,j,l] = exp(im*2*π*rand())*f([kx[i],ky[j],kz[l]],f_pars...)/sqrt(2.0*ω_val)
            g_t[i,j,l] = im*sqrt(2.0*ω_val)*g[i,j,l]
            ω_val_p = ω([kx_p[i];ky_p[j];kz_p[l]],m2)
            g_p[i,j,l] = f([kx_p[i],ky_p[j],kz_p[l]],f_pars)*exp(im*2*π*rand())/sqrt(2.0*ω_val_p)
            end
        end
    end
    return g, g_t, g_p
end

function policut(x,x0,r0,p)
    r2 = (x-x0)'*(x-x0)
    if (r0 - r2) >= 0.0
        return (r0 - r2)^p/r0^p
    else
        return 0.0
    end
end

function poli_step(x,x0,x1,p)
    if x < 0.0 || x0 < 0.0 || x1 < x0
        error("the variable is suposed to be positive")
    end
    if x0 < x && x < x0 + (x1-x0)/2.0
        return 1.0 - ((x-x0)^p/abs((x1-x0)/2.0)^p)/2.0 #+ (x1-x)^p/(x1-x0)^p
    elseif x >= x0 + (x1-x0)/2.0 && x < x1
        return +(x1-x)^p/abs((x1-x0)/2.0)^p/2.0
    elseif x <= x0
        return 1.0
    else
        return 0.0
    end
end

function step_cut(g,Box,J,p,percent,δ)
    x,y,z = get_coords(Box,J)
    g_cut = copy(g)
    x0 = (Box[2:2:end].+Box[1:2:end])./2.0
    L = (Box[2:2:end].-Box[1:2:end])
    r0 = minimum(abs.(L))*percent
    for i in 1:J[1]
        for j in 1:J[2]
            for k in 1:J[3]
                r = sqrt((x[i]-x0[1])^2+(y[j]-x0[2])^2+(z[k]-x0[3])^2)
                g_cut[i,j,k] = g[i,j,k] * poli_step(r,r0*(1.0-δ),r0,p)
            end
        end
    end
    return g_cut
end

function polinomial_cut(g,Box,J,p,percent)
    dx = differentials(Box,J)
    x,y,z = get_coords(Box,J)
    g_cut = copy(g)
    x0 = (Box[2:2:end].+Box[1:2:end])./2.0
    L = (Box[2:2:end].-Box[1:2:end])
    r0 = minimum(abs.(L))*percent
    r0 = r0^2
    for i in 1:J[1]
        for j in 1:J[2]
            for k in 1:J[3]
                g_cut[i,j,k] = g[i,j,k] * policut([x[i];y[j];z[k]],x0,r0,p)
            end
        end
    end
    return g_cut
end

function get_coords(Box,J)
    dx = differentials(Box,J)
    x = [Box_x[1] + (i-1)*dx[1] for i in 1:J[1]]
    y = [Box_x[3] + (i-1)*dx[2] for i in 1:J[2]]
    z = [Box_x[5] + (i-1)*dx[3] for i in 1:J[3]]
    return x, y ,z
end
