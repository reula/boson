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


function get_rho_old(ϕ,ϕ_t,Box_x,J;accuracy_order=2)

    Dx = periodic_derivative_operator(derivative_order=1, accuracy_order=2, xmin=Box_x[1], xmax=Box_x[2], N=J[1])
    Dy = periodic_derivative_operator(derivative_order=1, accuracy_order=2, xmin=Box_x[3], xmax=Box_x[4], N=J[2])
    Dz = periodic_derivative_operator(derivative_order=1, accuracy_order=2, xmin=Box_x[3], xmax=Box_x[4], N=J[3])

    Dϕ = Array{ComplexF64}(undef, 3, J...)
    ρ = Array{Float64}(undef,J...)

    for i in 1:J[1]
        for j in 1:J[2]
            Dϕ[3,i,j,:] = Dz*ϕ[i,j,:]
        end
    end
    for j in 1:J[2]
        for k in 1:J[3]
            Dϕ[1,:,j,k] = Dx*ϕ[:,j,k]
        end
    end
    for i in 1:J[1]
        for k in 1:J[3]
            Dϕ[2,i,:,k] = Dy*ϕ[i,:,k]
        end
    end

    ρ = zeros(J...)

    for i in 1:J[1]
        for j in 1:J[2]
            for k in 1:J[3]
                for d in 1:3
                    ρ[i,j,k] = ρ[i,j,k] + real(Dϕ[d,i,j,k]*conj(Dϕ[d,i,j,k]))
                end
                ρ[i,j,k] = ρ[i,j,k] + real(ϕ_t[i,j,k]*conj(ϕ_t[i,j,k]))
            end
        end
    end
    return 2.0*ρ[:,:,:]
end


function get_rho(ϕ,ϕ_t,Box_x,J;accuracy_order=2)

    Dx = periodic_derivative_operator(derivative_order=1, accuracy_order=2, xmin=Box_x[1], xmax=Box_x[2], N=J[1])
    Dy = periodic_derivative_operator(derivative_order=1, accuracy_order=2, xmin=Box_x[3], xmax=Box_x[4], N=J[2])
    Dz = periodic_derivative_operator(derivative_order=1, accuracy_order=2, xmin=Box_x[3], xmax=Box_x[4], N=J[3])

    Dϕ = Array{ComplexF64}(undef, 3, J...)
    ρ = Array{Float64}(undef,J...)

    for i in 1:J[1]
        for j in 1:J[2]
            Dϕ[3,i,j,:] = Dz*ϕ[i,j,:]
        end
    end
    for j in 1:J[2]
        for k in 1:J[3]
            Dϕ[1,:,j,k] = Dx*ϕ[:,j,k]
        end
    end
    for i in 1:J[1]
        for k in 1:J[3]
            Dϕ[2,i,:,k] = Dy*ϕ[i,:,k]
        end
    end

    π2 = zeros(J...)
    ∇2 = zeros(J...)
    V = zeros(J...)

    for i in 1:J[1]
        for j in 1:J[2]
            for k in 1:J[3]
                for d in 1:3
                    ∇2[i,j,k] = ∇2[i,j,k] + real(Dϕ[d,i,j,k]*conj(Dϕ[d,i,j,k]))
                end
                π2[i,j,k] = real(ϕ_t[i,j,k]*conj(ϕ_t[i,j,k]))
                V[i,j,k] = real(ϕ[i,j,k]*conj(ϕ[i,j,k]))
            end
        end
    end
    return 2.0*π2[:,:,:], 2.0*∇2[:,:,:], 2.0*V[:,:,:]
end



function load_data_full_h5(data_name,Box,J,ϕ,ϕ_t,N_fields)
    x, y, z = get_coords(Box,J)
    π2, ∇2, V = get_rho(ϕ,ϕ_t,Box,J);
    h5file = h5open("$(data_name).h5", "w") do file
    write(file, "coord0", x)
    write(file, "coord1", y)
    write(file, "coord2", z)
    write(file, "nvars", [N_fields])
    write(file, "var0", real.(ϕ))
    write(file, "var1", imag.(ϕ))
    write(file, "var2", real.(ϕ_t))
    write(file, "var3", imag.(ϕ_t))
    write(file, "var4", π2)
    write(file, "var5", ∇2)
    write(file, "var6", V)
    end
end



"""
function to relax an initial data
"""
function F(u,t,p)
    x,y,z,dxu,dyu,dzu,Dx,Dy,Dz,D2x,D2y,D2z,d2,dx,ρ,J,par = p 
    a, b, τ = par 
    d2 .= 0.0
    for i in 1:J[1]
        for j in 1:J[2]
            dzu[1,i,j,:] = Dz*u[1,i,j,:]
            d2[i,j,:] = D2z*u[1,i,j,:]
        end
    end
    for i in 1:J[1]
        for k in 1:J[3]
            dyu[1,i,:,k] = Dy*u[1,i,:,k]
            d2[i,:,k] += D2y*u[1,i,:,k] 
        end
    end
    for j in 1:J[2]
        for k in 1:J[3]
            dxu[1,:,j,k] = Dx*u[1,:,j,k]
            d2[:,j,k] += D2x*u[1,:,j,k] 
        end
    end

    du[1,:,:,:] .= u[2,:,:,:]
    @. du[2,:,:,:] .= d2[:,:,:]  - τ * u[2,:,:,:] + ρ[:,:,:]*u[1,:,:,:]^5

    #boundaries (fases)

    for k in (1,J[3])
        for i in 1:J[1]
            for j in 1:J[2]
                r = sqrt(x[i]^2 + y[j]^2 + z[k]^2)
                u_r = (x[i]*dxu[1,i,j,k] + y[j]*dyu[1,i,j,k] + z[k]*dzu[1,i,j,k])/r
                du[1,i,j,k] += -(u_r + (u[1,i,j,k] - 1.0)/r)
                du[2,i,j,k] += - (u[2,i,j,k] + a*(u_r  + b*(u[1,i,j,k] .- 1.0)/r))/dx[3]*2.0
            end
        end
    end
    for j in (1,J[2])
        for i in 1:J[1]
            for k in 1:J[3]
                r = sqrt(x[i]^2 + y[j]^2 + z[k]^2)
                u_r = (x[i]*dxu[1,i,j,k] + y[j]*dyu[1,i,j,k] + z[k]*dzu[1,i,j,k])/r
                du[1,i,j,k] += -(u_r + (u[1,i,j,k] - 1.0)/r)
                du[2,i,j,k] += - (u[2,i,j,k] + a*(u_r + b*(u[1,i,j,k] .- 1.0))/r)/dx[2]*2.0
            end
        end
    end
    for i in (1,J[1])
        for j in 1:J[2]
            for k in 1:J[3]
                r = sqrt(x[i]^2 + y[j]^2 + z[k]^2)
                u_r = (x[i]*dxu[1,i,j,k] + y[j]*dyu[1,i,j,k] + z[k]*dzu[1,i,j,k])/r
                du[1,i,j,k] += -(u_r + (u[1,i,j,k] - 1.0)/r)
                du[2,i,j,k] += - (u[2,i,j,k] + a*(u_r + b*(u[1,i,j,k] .- 1.0))/r)/dx[1]*2.0
            end
        end
    end
    
    #boundaries (edges)
    
if false
    for k in (1,J[3])
        for i in (1,J[1])
            for j in 2:J[2]-1
                r = sqrt(x[i]^2 + y[j]^2 + z[k]^2)
                du[2,i,j,k] += - (u[2,i,j,k] + a*(z[k]*dzu[i,j,k] + 0.0*x[i]*dxu[i,j,k] + 0.0*y[j]*dyu[i,j,k] - b*(u[1,i,j,k] .- 1.0)/r))/dx[3]/2.0
            end
        end
    end
    for j in (1,J[2])
        for i in (1,J[1])
            for k in 2:J[3]-1
                r = sqrt(x[i]^2 + y[j]^2 + z[k]^2)
                du[2,i,j,k] += - (u[2,i,j,k] + a*(0.0*z[k]*dzu[i,j,k] + x[i]*dxu[i,j,k] + y[j]*dyu[i,j,k] - b*(u[1,i,j,k] .- 1.0))/r)/dx[2]/2.0
            end
        end
    end
    for i in 2:J[1]-1
        for j in (1,J[2])
            for k in (1,J[3])
                r = sqrt(x[i]^2 + y[j]^2 + z[k]^2)
                du[2,i,j,k] += - (u[2,i,j,k] + a*(z[k]*dzu[i,j,k] + x[i]*dxu[i,j,k] + y[j]*dyu[i,j,k] - b*(u[1,i,j,k] .- 1.0))/r)/dx[1]/2.0
            end
        end
    end
    
    #boundaries (corners)

    for k in (1,J[3])
        for i in (1,J[1])
            for j in (1,J[2]) 
                r = sqrt(x[i]^2 + y[j]^2 + z[k]^2)
                du[2,i,j,k] += - (u[2,i,j,k] + a*(z[k]*dzu[i,j,k] + x[i]*dxu[i,j,k] + y[j]*dyu[i,j,k] - b*(u[1,i,j,k] .- 1.0)/r))/dx[3]/2.0
            end
        end
    end
end


    return du[:,:,:,:]
end


function FC(u,t,p)
    x,y,z,dxu,dyu,dzu,Dx,Dy,Dz,D2x,D2y,D2z,d2,dx,ρ,J,par = p 
    a, b, τ = par 
    n = #[0.0,0.0] 
    n = [1.0, 0.0] #asymtotic conditions
    d2 .= 0.0
    for i in 1:J[1]
        for j in 1:J[2]
            for d in 1:2
                dzu[d,i,j,:] = Dz*u[d,i,j,:]
            end
            d2[i,j,:] = D2z*u[1,i,j,:]
        end
    end
    for i in 1:J[1]
        for k in 1:J[3]
            for d in 1:2
                dyu[d,i,:,k] = Dy*u[d,i,:,k]
            end
            d2[i,:,k] += D2y*u[1,i,:,k] 
        end
    end
    for j in 1:J[2]
        for k in 1:J[3]
            for d in 1:2
                dxu[d,:,j,k] = Dx*u[d,:,j,k]
            end
            d2[:,j,k] += D2x*u[1,:,j,k] 
        end
    end

    du[1,:,:,:] .= u[2,:,:,:]
    @. du[2,:,:,:] .= d2[:,:,:] - τ * u[2,:,:,:] - ρ[:,:,:]#*u[1,:,:,:]^5

    #boundaries (fases)

    for k in (1,J[3])
        for i in 1:J[1]
            for j in 1:J[2]
                r = sqrt(x[i]^2 + y[j]^2 + z[k]^2)
                xdu = x[i]*dxu[:,i,j,k] + y[j]*dyu[:,i,j,k] + z[k]*dzu[:,i,j,k]
                du[:,i,j,k] = -a*(xdu + b*(u[:,i,j,k] .- n[:]))/r
            end
        end
    end
    for j in (1,J[2])
        for i in 1:J[1]
            for k in 1:J[3]
                r = sqrt(x[i]^2 + y[j]^2 + z[k]^2)
                xdu = x[i]*dxu[:,i,j,k] + y[j]*dyu[:,i,j,k] + z[k]*dzu[:,i,j,k]
                du[:,i,j,k] = -a*(xdu + b*(u[:,i,j,k] .- n[:]))/r
            end
        end
    end
    for i in (1,J[1])
        for j in 1:J[2]
            for k in 1:J[3]
                r = sqrt(x[i]^2 + y[j]^2 + z[k]^2)
                xdu = x[i]*dxu[:,i,j,k] + y[j]*dyu[:,i,j,k] + z[k]*dzu[:,i,j,k]
                du[:,i,j,k] = -a*(xdu + b*(u[:,i,j,k] .- n[:]))/r
            end
        end
    end
    
    #boundaries (edges)
    
    for k in (1,J[3])
        for i in (1,J[1])
            for j in 2:J[2]-1
                r = sqrt(x[i]^2 + y[j]^2 + z[k]^2)
                xdu = x[i]*dxu[:,i,j,k] + y[j]*dyu[:,i,j,k] + z[k]*dzu[:,i,j,k]
                du[:,i,j,k] = -a*(xdu + b*(u[:,i,j,k] .- n[:]))/r
            end
        end
    end
    for j in (1,J[2])
        for i in (1,J[1])
            for k in 2:J[3]-1
                r = sqrt(x[i]^2 + y[j]^2 + z[k]^2)
                xdu = x[i]*dxu[:,i,j,k] + y[j]*dyu[:,i,j,k] + z[k]*dzu[:,i,j,k]
                du[:,i,j,k] = -a*(xdu + b*(u[:,i,j,k] .- n[:]))/r
            end
        end
    end
    for i in 2:J[1]-1
        for j in (1,J[2])
            for k in (1,J[3])
                r = sqrt(x[i]^2 + y[j]^2 + z[k]^2)
                xdu = x[i]*dxu[:,i,j,k] + y[j]*dyu[:,i,j,k] + z[k]*dzu[:,i,j,k]
                du[:,i,j,k] = -a*(xdu + b*(u[:,i,j,k] .- n[:]))/r
            end
        end
    end
    
    #boundaries (corners)

    for k in (1,J[3])
        for i in (1,J[1])
            for j in (1,J[2]) 
                r = sqrt(x[i]^2 + y[j]^2 + z[k]^2)
                xdu = x[i]*dxu[:,i,j,k] + y[j]*dyu[:,i,j,k] + z[k]*dzu[:,i,j,k]
                du[:,i,j,k] = -a*(xdu + b*(u[:,i,j,k] .- n[:]))/r
            end
        end
    end


    return du[:,:,:,:]
end

function chichon(x,x0,Box,r0,p)
    d = x - x0
    r02 = r0^2
    r2 = d'*d
    if r2 < r0^2 
        return (r02 - r2)^p/r02^p
    else
        return 0.0
    end
end

function carlos(x,x0,Box,r0,A0)
    d = x - x0
    r02 = r0^2
    r = sqrt(d'*d)
    return A0*(-6 + 4*r^2/r02)/r02*exp(-r^2/r02)
end

function carlossol(x,x0,Box,r0,A0)
    d = x - x0
    r02 = r0^2
    r = sqrt(d'*d)
    return 1 + A0*exp(-r^2/r02)
end



function get_source(f, par)
    @show x,y,z,x0,Box_x, r0, p, J = par
    m = zeros(J...)
    for i in 1:J[1]
        for j in 1:J[2]
            for k in 1:J[3]
                m[i,j,k] = f([x[i],y[j],z[k]],x0,Box_x,r0,p)
            end
        end
    end
    return m[:,:,:]
end

