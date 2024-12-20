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



function capa(k,c_pars)
    k0, δ = c_pars
    if abs(k'*k - k0^2) < δ^2
        return (1.0 + 0.0im)/(2π*(k0*2*δ)) # normalized so that the integral is 1
    else
        return 0.0 + 0.0im
    end
end

function bola(k,b_pars)
    k0, δ  = b_pars
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
            g[i,j,l] = exp(im*2*π*rand())*f([kx[i],ky[j],kz[l]],f_pars)/sqrt(2.0*ω_val)
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

"""
add some extra grid point
"""
function get_coords_large(Box,J)
    dx = differentials(Box,J)
    x = [Box_x[1] + (i-1)*dx[1] for i in 1:J[1]+1]
    y = [Box_x[3] + (i-1)*dx[2] for i in 1:J[2]+1]
    z = [Box_x[5] + (i-1)*dx[3] for i in 1:J[3]+1]
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
    return π2[:,:,:], ∇2[:,:,:], V[:,:,:]
end



function load_data_full_h5(data_name,Box,J,ϕ,ϕ_t,N_fields)
    x, y, z = get_coords_large(Box,J)
    π2, ∇2, V = get_rho(ϕ,ϕ_t,Box,J);
    h5file = h5open("$(data_name).h5", "w") do file
    write(file, "coord0", x)
    write(file, "coord1", y)
    write(file, "coord2", z)
    write(file, "nvars", [N_fields])
    write(file, "var0", add_last_points(real.(ϕ)))
    write(file, "var1", add_last_points(imag.(ϕ)))
    write(file, "var2", add_last_points(real.(ϕ_t)))
    write(file, "var3", add_last_points(imag.(ϕ_t)))
    write(file, "var4", add_last_points(π2))
    write(file, "var5", add_last_points(∇2))
    write(file, "var6", add_last_points(V))
    end
end

"""
Add one face of periodic arrays to make them symmetric, that is we add a row and a column with the first values
"""
function add_last_points(m)
    J = size(m)
    J_large = J .+ (1,1,1)
    m_l = zeros(J_large)
    m_l[1:J[1],1:J[2],1:J[3]] = m
    m_l[J[1]+1,1:J[2],1:J[3]] = m[1,:,:]
    m_l[1:J[1],J[2]+1,1:J[3]] = m[:,1,:]
    m_l[1:J[1],1:J[2],J[3]+1] = m[:,:,1]
    m_l[1:J[1],J[2]+1,J[3]+1] = m[:,1,1]
    m_l[J[1]+1,J[2]+1,1:J[3]] = m[1,1,:]
    m_l[J[1]+1,1:J[2],J[3]+1] = m[1,:,1]
    m_l[J[1]+1,J[2]+1,J[3]+1] = m[1,1,1]
    return m_l[:,:,:]
end



"""
function to relax an initial data
version with Penalties
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
    @. du[2,:,:,:] .= d2[:,:,:]  - τ * u[2,:,:,:] - ρ[:,:,:]#*u[1,:,:,:]^5

    #boundaries (fases)

    for k in (1,J[3])
        for i in 2:J[1]-1
            for j in 2:J[2]-1
                r = sqrt(x[i]^2 + y[j]^2 + z[k]^2)
                u_r = (x[i]*dxu[1,i,j,k] + y[j]*dyu[1,i,j,k] + z[k]*dzu[1,i,j,k])/r
                du[1,i,j,k] += -(u_r + (u[1,i,j,k] - 1.0)/r)
                du[2,i,j,k] += - (u[2,i,j,k] + a*(u_r  + b*(u[1,i,j,k] .- 1.0)/r))/dx[3]*2.0
            end
        end
    end
    for j in (1,J[2])
        for i in 2:J[1]-1
            for k in 2:J[3]-1
                r = sqrt(x[i]^2 + y[j]^2 + z[k]^2)
                u_r = (x[i]*dxu[1,i,j,k] + y[j]*dyu[1,i,j,k] + z[k]*dzu[1,i,j,k])/r
                du[1,i,j,k] += -(u_r + (u[1,i,j,k] - 1.0)/r)
                du[2,i,j,k] += - (u[2,i,j,k] + a*(u_r + b*(u[1,i,j,k] .- 1.0))/r)/dx[2]*2.0
            end
        end
    end
    for i in (1,J[1])
        for j in 2:J[2]-1
            for k in 2:J[3]-1
                r = sqrt(x[i]^2 + y[j]^2 + z[k]^2)
                u_r = (x[i]*dxu[1,i,j,k] + y[j]*dyu[1,i,j,k] + z[k]*dzu[1,i,j,k])/r
                du[1,i,j,k] += -(u_r + (u[1,i,j,k] - 1.0)/r)
                du[2,i,j,k] += - (u[2,i,j,k] + a*(u_r + b*(u[1,i,j,k] .- 1.0))/r)/dx[1]*2.0
            end
        end
    end
    
    #boundaries (edges)
    
if true
    for k in (1,J[3])
        for i in (1,J[1])
            for j in 2:J[2]-1
                r = sqrt(x[i]^2 + y[j]^2 + z[k]^2)
                u_r = (x[i]*dxu[1,i,j,k] + y[j]*dyu[1,i,j,k] + z[k]*dzu[1,i,j,k])/r
                du[1,i,j,k] += -(u_r + (u[1,i,j,k] - 1.0)/r)
                du[2,i,j,k] += - (u[2,i,j,k] + a*(u_r + b*(u[1,i,j,k] .- 1.0))/r)/dx[1]*2.0
            end
        end
    end
    for j in (1,J[2])
        for i in (1,J[1])
            for k in 2:J[3]-1
                r = sqrt(x[i]^2 + y[j]^2 + z[k]^2)
                u_r = (x[i]*dxu[1,i,j,k] + y[j]*dyu[1,i,j,k] + z[k]*dzu[1,i,j,k])/r
                du[1,i,j,k] += -(u_r + (u[1,i,j,k] - 1.0)/r)
                du[2,i,j,k] += - (u[2,i,j,k] + a*(u_r + b*(u[1,i,j,k] .- 1.0))/r)/dx[1]*2.0
            end
        end
    end
    for i in 2:J[1]-1
        for j in (1,J[2])
            for k in (1,J[3])
                r = sqrt(x[i]^2 + y[j]^2 + z[k]^2)
                u_r = (x[i]*dxu[1,i,j,k] + y[j]*dyu[1,i,j,k] + z[k]*dzu[1,i,j,k])/r
                du[1,i,j,k] += -(u_r + (u[1,i,j,k] - 1.0)/r)
                du[2,i,j,k] += - (u[2,i,j,k] + a*(u_r + b*(u[1,i,j,k] .- 1.0))/r)/dx[1]*2.0
            end
        end
    end
    
    #boundaries (corners)

    for k in (1,J[3])
        for i in (1,J[1])
            for j in (1,J[2]) 
                r = sqrt(x[i]^2 + y[j]^2 + z[k]^2)
                u_r = (x[i]*dxu[1,i,j,k] + y[j]*dyu[1,i,j,k] + z[k]*dzu[1,i,j,k])/r
                du[1,i,j,k] += -(u_r + (u[1,i,j,k] - 1.0)/r)
                du[2,i,j,k] += - (u[2,i,j,k] + a*(u_r + b*(u[1,i,j,k] .- 1.0))/r)/dx[1]*2.0
            end
        end
    end
end


    return du[:,:,:,:]
end
"""
function to relax initial data, version with boudary conditions imposed at boundary explicitly.
"""

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

"""
function to relax initial data, version with boudary conditions imposed at boundary explicitly.
"""
function FCR(u,t,p)
    coords,boundary_derivs,derivs,d2,dx,ρ,J,par = p 
    x,y,z = coords
    dxu_x,dyu_x,dzu_x,dxu_y,dyu_y,dzu_y,dxu_z,dyu_z,dzu_z = boundary_derivs 
    Dx,Dy,Dz,D2x,D2y,D2z = derivs
    a, b, τ = par 
    n = #[0.0,0.0] 
    n = [1.0, 0.0] #asymtotic conditions
    d2 .= 0.0
#face z (i,j)
    for i in 1:J[1]
        for j in 1:J[2]
            for d in 1:2
                dzu_z[1,d,i,j] = derivative_left(D2z, u[d,i,j,:], Val{1}())
                dzu_z[2,d,i,j] = derivative_right(D2z, u[d,i,j,:], Val{1}())
            end
            d2[i,j,:] = D2z*u[1,i,j,:]
        end
        for d in 1:2
            dyu_z[1,d,i,:] = Dy*u[d,i,:,1]
            dyu_z[2,d,i,:] = Dy*u[d,i,:,J[3]]
        end
    end
    for j in 1:J[2]
        for d in 1:2
            dxu_z[1,d,:,j] = Dx*u[d,:,j,1]
            dxu_z[2,d,:,j] = Dx*u[d,:,j,J[3]]
        end
    end
#face y (i,k)
    for i in 1:J[1]
        for k in 1:J[3]
            for d in 1:2
                dzu_y[1,d,i,k] = derivative_left(D2y,  u[d,i,:,k], Val{1}())
                dzu_y[2,d,i,k] = derivative_right(D2y, u[d,i,:,k], Val{1}())
            end
            d2[i,:,k] += D2y*u[1,i,:,k]
        end
        for d in 1:2
            dzu_y[1,d,i,:] = Dz*u[d,i,1,:]
            dzu_y[2,d,i,:] = Dz*u[d,i,J[2],:]
        end
    end
    for k in 1:J[3]
        for d in 1:2
            dxu_y[1,d,:,k] = Dx*u[d,:,1,k]
            dxu_y[2,d,:,k] = Dx*u[d,:,J[2],k]
        end
    end
#face x (j,k)
for j in 1:J[2]
    for k in 1:J[3]
        for d in 1:2
            dzu_x[1,d,j,k] = derivative_left(D2x,  u[d,:,j,k], Val{1}())
            dzu_x[2,d,j,k] = derivative_right(D2x, u[d,:,j,k], Val{1}())
        end
        d2[:,j,k] += D2x*u[1,:,j,k]
    end
    for d in 1:2
        dzu_x[1,d,j,:] = Dz*u[d,1,j,:]
        dzu_x[2,d,j,:] = Dz*u[d,J[1],j,:]
    end
end
for k in 1:J[3]
    for d in 1:2
        dyu_x[1,d,:,k] = Dx*u[d,1,:,k]
        dyu_x[2,d,:,k] = Dx*u[d,J[2],:,k]
    end
end

    du[1,:,:,:] .= u[2,:,:,:]
    @. du[2,:,:,:] .= d2[:,:,:] - τ * u[2,:,:,:] - ρ[:,:,:]#*u[1,:,:,:]^5

    #boundaries (fases)

    for k in (1,J[3])
        for i in 1:J[1]
            for j in 1:J[2]
                r = sqrt(x[i]^2 + y[j]^2 + z[k]^2)
                xdu = x[i]*dxu_z[fase(k),:,i,j] + y[j]*dyu_z[fase(k),:,i,j] + z[k]*dzu_z[fase(k),:,i,j]
                du[:,i,j,k] = -a*(xdu + b*(u[:,i,j,k] .- n[:]))/r
            end
        end
    end
    for j in (1,J[2])
        for i in 1:J[1]
            for k in 1:J[3]
                r = sqrt(x[i]^2 + y[j]^2 + z[k]^2)
                xdu = x[i]*dxu_y[fase(j),:,i,k] + y[j]*dyu_y[fase(j),:,i,k] + z[k]*dzu_y[fase(j),:,i,k]
                du[:,i,j,k] = -a*(xdu + b*(u[:,i,j,k] .- n[:]))/r
            end
        end
    end
    for i in (1,J[1])
        for j in 1:J[2]
            for k in 1:J[3]
                r = sqrt(x[i]^2 + y[j]^2 + z[k]^2)
                xdu = x[i]*dxu_x[fase(i),:,j,k] + y[j]*dyu_x[fase(i),:,j,k] + z[k]*dzu_x[fase(i),:,j,k]
                du[:,i,j,k] = -a*(xdu + b*(u[:,i,j,k] .- n[:]))/r
            end
        end
    end
    
    #boundaries (edges)
    
    for k in (1,J[3])
        for i in (1,J[1])
            for j in 2:J[2]-1
                r = sqrt(x[i]^2 + y[j]^2 + z[k]^2)
                xdu = x[i]*dxu_z[fase(k),:,i,j] + y[j]*dyu_z[fase(k),:,i,j] + z[k]*dzu_z[fase(k),:,i,j]
                du[:,i,j,k] = -a*(xdu + b*(u[:,i,j,k] .- n[:]))/r
            end
        end
    end
    for j in (1,J[2])
        for i in (1,J[1])
            for k in 2:J[3]-1
                r = sqrt(x[i]^2 + y[j]^2 + z[k]^2)
                xdu = x[i]*dxu_y[fase(j),:,i,k] + y[j]*dyu_y[fase(j),:,i,k] + z[k]*dzu_y[fase(k),:,i,j]
                du[:,i,j,k] = -a*(xdu + b*(u[:,i,j,k] .- n[:]))/r
            end
        end
    end
    for i in 2:J[1]-1
        for j in (1,J[2])
            for k in (1,J[3])
                r = sqrt(x[i]^2 + y[j]^2 + z[k]^2)
                xdu = x[i]*dxu_x[fase(i),:,j,k] + y[j]*dyu_x[fase(i),:,j,k] + z[k]*dzu_x[fase(i),:,j,k]
                du[:,i,j,k] = -a*(xdu + b*(u[:,i,j,k] .- n[:]))/r
            end
        end
    end
    
    #boundaries (corners)

    for k in (1,J[3])
        for i in (1,J[1])
            for j in (1,J[2]) 
                r = sqrt(x[i]^2 + y[j]^2 + z[k]^2)
                xdu = x[i]*dxu_x[fase(i),:,j,k] + y[j]*dyu_x[fase(i),:,j,k] + z[k]*dzu_x[fase(i),:,j,k]
                du[:,i,j,k] = -a*(xdu + b*(u[:,i,j,k] .- n[:]))/r
            end
        end
    end


    return du[:,:,:,:]
end

"""
function to relax initial data and lapse, version with boudary conditions imposed at boundary explicitly.
"""

function FCR_Full(u,t,p)
    coords,boundary_derivs,derivs,d2,dx,ρ,J,n_fields,par = p 
    x,y,z = coords
    dxu_x,dyu_x,dzu_x,dxu_y,dyu_y,dzu_y,dxu_z,dyu_z,dzu_z = boundary_derivs 
    Dx,Dy,Dz,D2x,D2y,D2z = derivs
    π2,∇2,V = ρ
    a, b, τ = par 
    n = #[0.0,0.0] 
    n = [1.0, 0.0,1.0,0.0] #asymtotic conditions
    d2 .= 0.0
#face z (i,j)
Threads.@threads for i in 1:J[1]
        for j in 1:J[2]
            for d in 1:n_fields
                dzu_z[1,d,i,j] = derivative_left(D2z, u[d,i,j,:], Val{1}())
                dzu_z[2,d,i,j] = derivative_right(D2z, u[d,i,j,:], Val{1}())
            end
            d2[1,i,j,:] = D2z*u[1,i,j,:]
            d2[2,i,j,:] = D2z*u[3,i,j,:]
        end
        for d in 1:n_fields
            dyu_z[1,d,i,:] = Dy*u[d,i,:,1]
            dyu_z[2,d,i,:] = Dy*u[d,i,:,J[3]]
        end
    end
    for j in 1:J[2]
        for d in 1:n_fields
            dxu_z[1,d,:,j] = Dx*u[d,:,j,1]
            dxu_z[2,d,:,j] = Dx*u[d,:,j,J[3]]
        end
    end
#face y (i,k)
Threads.@threads    for i in 1:J[1]
        for k in 1:J[3]
            for d in 1:n_fields
                dyu_y[1,d,i,k] = derivative_left(D2y,  u[d,i,:,k], Val{1}())
                dyu_y[2,d,i,k] = derivative_right(D2y, u[d,i,:,k], Val{1}())
            end
            d2[1,i,:,k] += D2y*u[1,i,:,k]
            d2[2,i,:,k] += D2y*u[3,i,:,k]
        end
        for d in 1:n_fields
            dzu_y[1,d,i,:] = Dz*u[d,i,1,:]
            dzu_y[2,d,i,:] = Dz*u[d,i,J[2],:]
        end
    end
    for k in 1:J[3]
        for d in 1:n_fields
            dxu_y[1,d,:,k] = Dx*u[d,:,1,k]
            dxu_y[2,d,:,k] = Dx*u[d,:,J[2],k]
        end
    end
#face x (j,k)
Threads.@threads for j in 1:J[2]
    for k in 1:J[3]
        for d in 1:n_fields
            dxu_x[1,d,j,k] = derivative_left(D2x,  u[d,:,j,k], Val{1}())
            dxu_x[2,d,j,k] = derivative_right(D2x, u[d,:,j,k], Val{1}())
        end
        d2[1,:,j,k] += D2x*u[1,:,j,k]
        d2[2,:,j,k] += D2x*u[3,:,j,k]
    end
    for d in 1:n_fields
        dzu_x[1,d,j,:] = Dz*u[d,1,j,:]
        dzu_x[2,d,j,:] = Dz*u[d,J[1],j,:]
    end
end
for k in 1:J[3]
    for d in 1:n_fields
        dyu_x[1,d,:,k] = Dy*u[d,1,:,k]
        dyu_x[2,d,:,k] = Dy*u[d,J[1],:,k]
    end
end

    du[1,:,:,:] .= u[2,:,:,:]
    du[3,:,:,:] .= u[4,:,:,:]
    
    @. du[2,:,:,:] .= d2[1,:,:,:] - τ * u[2,:,:,:] +  2π*((π2[:,:,:] + V[:,:,:])*u[1,:,:,:]^5 + ∇2[:,:,:]*u[1,:,:,:])
    @. du[4,:,:,:] .= d2[2,:,:,:] - τ * u[4,:,:,:] -  2π*u[3,:,:,:]*((u[1,:,:,:]^4)*(7*π2[:,:,:] - 5*V[:,:,:]) - ∇2[:,:,:])

    #boundaries (fases)

    for k in (1,J[3])
        for i in 1:J[1]
            for j in 1:J[2]
                r = sqrt(x[i]^2 + y[j]^2 + z[k]^2)
                xdu = x[i]*dxu_z[fase(k),:,i,j] + y[j]*dyu_z[fase(k),:,i,j] + z[k]*dzu_z[fase(k),:,i,j]
                du[:,i,j,k] = -a*(xdu + b*(u[:,i,j,k] .- n[:]))/r
            end
        end
    end
    for j in (1,J[2])
        for i in 1:J[1]
            for k in 1:J[3]
                r = sqrt(x[i]^2 + y[j]^2 + z[k]^2)
                xdu = x[i]*dxu_y[fase(j),:,i,k] + y[j]*dyu_y[fase(j),:,i,k] + z[k]*dzu_y[fase(j),:,i,k]
                du[:,i,j,k] = -a*(xdu + b*(u[:,i,j,k] .- n[:]))/r
            end
        end
    end
    for i in (1,J[1])
        for j in 1:J[2]
            for k in 1:J[3]
                r = sqrt(x[i]^2 + y[j]^2 + z[k]^2)
                xdu = x[i]*dxu_x[fase(i),:,j,k] + y[j]*dyu_x[fase(i),:,j,k] + z[k]*dzu_x[fase(i),:,j,k]
                du[:,i,j,k] = -a*(xdu + b*(u[:,i,j,k] .- n[:]))/r
            end
        end
    end
    
    #boundaries (edges)
    
    for k in (1,J[3])
        for i in (1,J[1])
            for j in 2:J[2]-1
                r = sqrt(x[i]^2 + y[j]^2 + z[k]^2)
                xdu = x[i]*dxu_z[fase(k),:,i,j] + y[j]*dyu_z[fase(k),:,i,j] + z[k]*dzu_z[fase(k),:,i,j]
                du[:,i,j,k] = -a*(xdu + b*(u[:,i,j,k] .- n[:]))/r
            end
        end
    end
    for j in (1,J[2])
        for i in (1,J[1])
            for k in 2:J[3]-1
                r = sqrt(x[i]^2 + y[j]^2 + z[k]^2)
                xdu = x[i]*dxu_y[fase(j),:,i,k] + y[j]*dyu_y[fase(j),:,i,k] + z[k]*dzu_y[fase(k),:,i,j]
                du[:,i,j,k] = -a*(xdu + b*(u[:,i,j,k] .- n[:]))/r
            end
        end
    end
    for i in 2:J[1]-1
        for j in (1,J[2])
            for k in (1,J[3])
                r = sqrt(x[i]^2 + y[j]^2 + z[k]^2)
                xdu = x[i]*dxu_x[fase(i),:,j,k] + y[j]*dyu_x[fase(i),:,j,k] + z[k]*dzu_x[fase(i),:,j,k]
                du[:,i,j,k] = -a*(xdu + b*(u[:,i,j,k] .- n[:]))/r
            end
        end
    end
    
    #boundaries (corners)

    for k in (1,J[3])
        for i in (1,J[1])
            for j in (1,J[2]) 
                r = sqrt(x[i]^2 + y[j]^2 + z[k]^2)
                xdu = x[i]*dxu_x[fase(i),:,j,k] + y[j]*dyu_x[fase(i),:,j,k] + z[k]*dzu_x[fase(i),:,j,k]
                du[:,i,j,k] = -a*(xdu + b*(u[:,i,j,k] .- n[:]))/r
            end
        end
    end


    return du[:,:,:,:]
end

"""
function to relax initial data and lapse, version with boudary conditions imposed with Penalties.
"""
function F_full(u,t,p)
    coords,boundary_derivs,derivs,d2,dx,ρ,J,n_fields,par = p 
    x,y,z = coords
    dxu_x,dyu_x,dzu_x,dxu_y,dyu_y,dzu_y,dxu_z,dyu_z,dzu_z = boundary_derivs 
    Dx,Dy,Dz,D2x,D2y,D2z = derivs
    π2,∇2,V = ρ
    a, b, τ = par 
    n = #[0.0,0.0] 
    n = [1.0, 0.0,1.0,0.0] #asymtotic conditions
    d2 .= 0.0
#face z (i,j)
Threads.@threads for i in 1:J[1]
        for j in 1:J[2]
            for d in 1:n_fields
                dzu_z[1,d,i,j] = derivative_left(D2z, u[d,i,j,:], Val{1}())
                dzu_z[2,d,i,j] = derivative_right(D2z, u[d,i,j,:], Val{1}())
            end
            d2[1,i,j,:] = D2z*u[1,i,j,:]
            d2[2,i,j,:] = D2z*u[3,i,j,:]
        end
        for d in 1:n_fields
            dyu_z[1,d,i,:] = Dy*u[d,i,:,1]
            dyu_z[2,d,i,:] = Dy*u[d,i,:,J[3]]
        end
    end
    for j in 1:J[2]
        for d in 1:n_fields
            dxu_z[1,d,:,j] = Dx*u[d,:,j,1]
            dxu_z[2,d,:,j] = Dx*u[d,:,j,J[3]]
        end
    end
#face y (i,k)
Threads.@threads    for i in 1:J[1]
        for k in 1:J[3]
            for d in 1:n_fields
                dzu_y[1,d,i,k] = derivative_left(D2y,  u[d,i,:,k], Val{1}())
                dzu_y[2,d,i,k] = derivative_right(D2y, u[d,i,:,k], Val{1}())
            end
            d2[1,i,:,k] += D2y*u[1,i,:,k]
            d2[2,i,:,k] += D2y*u[3,i,:,k]
        end
        for d in 1:n_fields
            dzu_y[1,d,i,:] = Dz*u[d,i,1,:]
            dzu_y[2,d,i,:] = Dz*u[d,i,J[2],:]
        end
    end
    for k in 1:J[3]
        for d in 1:n_fields
            dxu_y[1,d,:,k] = Dx*u[d,:,1,k]
            dxu_y[2,d,:,k] = Dx*u[d,:,J[2],k]
        end
    end
#face x (j,k)
Threads.@threads for j in 1:J[2]
    for k in 1:J[3]
        for d in 1:n_fields
            dxu_x[1,d,j,k] = derivative_left(D2x,  u[d,:,j,k], Val{1}())
            dxu_x[2,d,j,k] = derivative_right(D2x, u[d,:,j,k], Val{1}())
        end
        d2[1,:,j,k] += D2x*u[1,:,j,k]
        d2[2,:,j,k] += D2x*u[3,:,j,k]
    end
    for d in 1:n_fields
        dzu_x[1,d,j,:] = Dz*u[d,1,j,:]
        dzu_x[2,d,j,:] = Dz*u[d,J[1],j,:]
    end
end
for k in 1:J[3]
    for d in 1:n_fields
        dyu_x[1,d,:,k] = Dy*u[d,1,:,k]
        dyu_x[2,d,:,k] = Dy*u[d,J[1],:,k]
    end
end

    du[1,:,:,:] .= u[2,:,:,:]
    du[3,:,:,:] .= u[4,:,:,:]
    
    @. du[2,:,:,:] .= d2[1,:,:,:] - τ * u[2,:,:,:] +  2π*((π2[:,:,:] + V[:,:,:])*u[1,:,:,:]^5 + ∇2[:,:,:]*u[1,:,:,:])
    @. du[4,:,:,:] .= d2[2,:,:,:] - τ * u[4,:,:,:] -  2π*u[3,:,:,:]*((u[1,:,:,:]^4)*(7*π2[:,:,:] - 5*V[:,:,:]) - ∇2[:,:,:])

    #boundaries (fases)

    for k in (1,J[3])
        for i in 1:J[1]
            for j in 1:J[2]
                r = sqrt(x[i]^2 + y[j]^2 + z[k]^2)
                xdu = x[i]*dxu_z[fase(k),:,i,j] + y[j]*dyu_z[fase(k),:,i,j] + z[k]*dzu_z[fase(k),:,i,j]
                du[1,i,j,k] += -(xdu[1] + (u[1,i,j,k] .- n[1]))/r
                du[2,i,j,k] += -(u[2,i,j,k] + (xdu[1] + (u[1,i,j,k] .- n[1]))/r)/dx[3]*2.0
                du[3,i,j,k] += -(xdu[3] + (u[3,i,j,k] .- n[3]))/r
                du[4,i,j,k] += -(u[4,i,j,k] + (xdu[3] + (u[3,i,j,k] .- n[3]))/r)/dx[3]*2.0
            end
        end
    end
    for j in (1,J[2])
        for i in 1:J[1]
            for k in 1:J[3]
                r = sqrt(x[i]^2 + y[j]^2 + z[k]^2)
                xdu = x[i]*dxu_y[fase(j),:,i,k] + y[j]*dyu_y[fase(j),:,i,k] + z[k]*dzu_y[fase(j),:,i,k]
                du[1,i,j,k] += -(xdu[1] + (u[1,i,j,k] .- n[1]))/r
                du[2,i,j,k] += -(u[2,i,j,k] + (xdu[1] + (u[1,i,j,k] .- n[1]))/r)/dx[2]*2.0
                du[3,i,j,k] += -(xdu[3] + (u[3,i,j,k] .- n[3]))/r
                du[4,i,j,k] += -(u[4,i,j,k] + (xdu[3] + (u[3,i,j,k] .- n[3]))/r)/dx[2]*2.0
            end
        end
    end
    for i in (1,J[1])
        for j in 1:J[2]
            for k in 1:J[3]
                r = sqrt(x[i]^2 + y[j]^2 + z[k]^2)
                xdu = x[i]*dxu_x[fase(i),:,j,k] + y[j]*dyu_x[fase(i),:,j,k] + z[k]*dzu_x[fase(i),:,j,k]
                du[1,i,j,k] += -(xdu[1] + (u[1,i,j,k] .- n[1]))/r
                du[2,i,j,k] += -(u[2,i,j,k] + (xdu[1] + (u[1,i,j,k] .- n[1]))/r)/dx[1]*2.0
                du[3,i,j,k] += -(xdu[3] + (u[3,i,j,k] .- n[3]))/r
                du[4,i,j,k] += -(u[4,i,j,k] + (xdu[3] + (u[3,i,j,k] .- n[3]))/r)/dx[1]*2.0
            end
        end
    end
    
    #boundaries (edges)
    
    for k in (1,J[3])
        for i in (1,J[1])
            for j in 2:J[2]-1
                r = sqrt(x[i]^2 + y[j]^2 + z[k]^2)
                xdu = x[i]*dxu_z[fase(k),:,i,j] + y[j]*dyu_z[fase(k),:,i,j] + z[k]*dzu_z[fase(k),:,i,j]
                du[1,i,j,k] += -(xdu[1] + (u[1,i,j,k] .- n[1]))/r
                du[2,i,j,k] += -(u[2,i,j,k] + (xdu[1] + (u[1,i,j,k] .- n[1]))/r)/dx[1]*2.0
                du[3,i,j,k] += -(xdu[3] + (u[3,i,j,k] .- n[3]))/r
                du[4,i,j,k] += -(u[4,i,j,k] + (xdu[3] + (u[3,i,j,k] .- n[3]))/r)/dx[1]*2.0
            end
        end
    end
    for j in (1,J[2])
        for i in (1,J[1])
            for k in 2:J[3]-1
                r = sqrt(x[i]^2 + y[j]^2 + z[k]^2)
                xdu = x[i]*dxu_y[fase(j),:,i,k] + y[j]*dyu_y[fase(j),:,i,k] + z[k]*dzu_y[fase(k),:,i,j]
                du[1,i,j,k] += -(xdu[1] + (u[1,i,j,k] .- n[1]))/r
                du[2,i,j,k] += -(u[2,i,j,k] + (xdu[1] + (u[1,i,j,k] .- n[1]))/r)/dx[1]*2.0
                du[3,i,j,k] += -(xdu[3] + (u[3,i,j,k] .- n[3]))/r
                du[4,i,j,k] += -(u[4,i,j,k] + (xdu[3] + (u[3,i,j,k] .- n[3]))/r)/dx[1]*2.0
            end
        end
    end
    for i in 2:J[1]-1
        for j in (1,J[2])
            for k in (1,J[3])
                r = sqrt(x[i]^2 + y[j]^2 + z[k]^2)
                xdu = x[i]*dxu_x[fase(i),:,j,k] + y[j]*dyu_x[fase(i),:,j,k] + z[k]*dzu_x[fase(i),:,j,k]
                du[1,i,j,k] += -(xdu[1] + (u[1,i,j,k] .- n[1]))/r
                du[2,i,j,k] += -(u[2,i,j,k] + (xdu[1] + (u[1,i,j,k] .- n[1]))/r)/dx[1]*2.0
                du[3,i,j,k] += -(xdu[3] + (u[3,i,j,k] .- n[3]))/r
                du[4,i,j,k] += -(u[4,i,j,k] + (xdu[3] + (u[3,i,j,k] .- n[3]))/r)/dx[1]*2.0
              end
        end
    end
    
    #boundaries (corners)

    for k in (1,J[3])
        for i in (1,J[1])
            for j in (1,J[2]) 
                r = sqrt(x[i]^2 + y[j]^2 + z[k]^2)
                xdu = x[i]*dxu_x[fase(i),:,j,k] + y[j]*dyu_x[fase(i),:,j,k] + z[k]*dzu_x[fase(i),:,j,k]
                du[1,i,j,k] += -(xdu[1] + (u[1,i,j,k] .- n[1]))/r
                du[2,i,j,k] += -(u[2,i,j,k] + (xdu[1] + (u[1,i,j,k] .- n[1]))/r)/dx[1]*2.0
                du[3,i,j,k] += -(xdu[3] + (u[3,i,j,k] .- n[3]))/r
                du[4,i,j,k] += -(u[4,i,j,k] + (xdu[3] + (u[3,i,j,k] .- n[3]))/r)/dx[1]*2.0
             end
        end
    end


    return du[:,:,:,:]
end


function fase(k) 
    if k == 1 
        return 1
    else
        return 2
    end
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
    x,y,z,x0,Box_x, r0, p, J = par
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

function embed_source(m,J,enlarge_factor)
    J_in = size(m)
    if enlarge_factor == 3
        if J == J_in .* 3
            m_l = zeros(J...)
            m_l[(J_in[1]+1):2*J_in[1],(J_in[2]+1):2*J_in[2],(J_in[3]+1):2*J_in[3]] = m[:,:,:]
        else
            error("J sizes don't match")
        end
    elseif enlarge_factor == 2
        if J == J_in .+ 2 .*((J_in .÷2) .+1)
            m_l = zeros(J...)
            m_l[J_in[1]÷2+2:J_in[1]÷2+1+J_in[1],J_in[2]÷2+2:J_in[2]÷2+1+J_in[2],J_in[3]÷2+2:J_in[3]÷2+1+J_in[3]] = m[:,:,:]
        else
            error("J sizes don't match J = $J, computed = $(J_in .+ 2 .*((J_in .÷2) .+1))")
        end
    else
        error("enlarge_factor = $(enlarge_factor) not implemented")
    end
    return m_l[:,:,:]
end

get_norm_time(v,m) = [norm(v[m,field,:,:,:]) for field in 1:lastindex(v[1,:,1,1,1])]/sqrt(prod(size(v[1,1,:,:,:])))
get_norms(u) = [norm(u[field,:,:,:]) for field in 1:lastindex(u[:,1,1,1])]/sqrt(prod(size(u[1,:,:,:])))


#=
@. model_phi(r,p) = (1 + p[1]/2/r)
@. model_v(r,pv) = (1 - pv[1]/2/r)

fit_phi = curve_fit(model_phi, x[200:J[1]], v[m,1,200:J[1],J[2] ÷ 2 , J[3]÷ 2], p0) 
=#

function get_coords(Box, J; periodic = true)
    D = length(J)
    d = 1
    x = [Box[2*d-1] + (Box[2d]-Box[2*d-1])/J[d]*(i-1) for i in 1:J[d]]
    if D == 1
        return x
    elseif D == 2
        d = 2
        y = [Box[2*d-1] + (Box[2d]-Box[2*d-1])/J[d]*(i-1) for i in 1:J[d]]
        return x, y
    elseif D == 3
        d = 2
        y = [Box[2*d-1] + (Box[2d]-Box[2*d-1])/J[d]*(i-1) for i in 1:J[d]]
        d = 3
        z = [Box[2*d-1] + (Box[2d]-Box[2*d-1])/J[d]*(i-1) for i in 1:J[d]]
        return x, y, z 
    else 
        error("not implemented for D > 3, D = $D")
    end
end