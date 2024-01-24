module BeamEstimation

using  NFFT, Unitful, AbstractFFTs, Requires, LinearAlgebra, ArrayTools
import UnitfulAngles: mas

export 	dirtybeam,
		maxpixsize,
		maxfov,
		beamellipse,
		mas


maxpixsize_(Bmax,λmin) = λmin/Bmax/2 
"""
    maxpixsize(Bmax,λ)

Compute the largest possible pixel size to prevent aliasing given :
-  `λ` the smallest wavelenght. If `λ` is an Array, it is computed on the minimum of `λ`
-  `Bmax` the largest baseline length. If `Bmax` is an Matrix, the second dimension encodes for the lateral dimension of the baseline
largest 

If `λ` and `Bmax` are unitless it is assume that they are in meters
"""
maxpixsize(Bmax::Real,λmin::Real) = maxpixsize_(Bmax,λmin) 
maxpixsize(Bmax::Unitful.Length,λmin::Unitful.Length)  = (maxpixsize_(Bmax,λmin)  |> NoUnits) *1u"rad" |> mas
maxpixsize(Bmax,λ::AbstractArray) = maxpixsize(Bmax,minimum(λ))
maxpixsize(baselines::AbstractMatrix,λ) = maxpixsize(maximum(sqrt.(sum(abs2,baselines,dims=2))),λ)


"""
    maxpixsize(baselinesU::AbstractVector,baselinesV::AbstractVector,λ)

	Compute the largest possible pixel size to prevent aliasing given :
	-  `λ` the smallest wavelenght. If `λmin` is an Array, it is computed on the minimum of `λmin`
	-  `baselinesU`  and `baselinesV` are the baseline legnth along each lateral dimensions  
"""
maxpixsize(baselinesU::AbstractVector,baselinesV::AbstractVector,λ) = maxpixsize(λ,maximum(sqrt.(abs2.(baselinesU) .+ abs2.(baselinesV))))

"""
    maxpixsize(uv) 
	Compute the largest possible pixel size to prevent aliasing  from the (uv) coordinates `uv`

"""
maxpixsize(uv::AbstractArray) = minimum((sum(abs2,uv,dims=1)).^(-1/2)) |> mas
maxpixsize(uv::AbstractArray{T,N}) where {N,T<:Unitful.Quantity} = minimum((sum(abs2,uv,dims=1)).^(-1/2)) / 2 |> mas

maxfov_(D,λmax) = λmax/D
"""
    maxfov(D,λ)

Computed the maximum field of view of the interferometer given the wavelenghts `λ` and the telescope diameter `D`.
If `λ` and `D` are unitless it is assume that they are in meters.
"""
maxfov(D::Real,λmax::Real) = maxfov_(D,λmax) 
maxfov(D::Unitful.Length,λmax::Unitful.Length) = (maxfov_(D,λmax)  |> NoUnits) *1u"rad" |> mas
maxfov(D,λ::AbstractArray) = maxfov(D,minimum(λ))

uvplane_(Baselines, λ) = Baselines ./ reshape(λ, 1,1,:)
function uvplane_(::Type{T},Baselines, λ) where T<:Unitful.Length
    uvplane_(Baselines, λ)  .|> u"rad^-1"
end
function uvplane_(::Type{T},Baselines, λ) where T<:Real
    uvplane_(Baselines, λ)  
end

"""
    uvplane(Baselines, λ)

Compute the uv-plane from the interferometer array baselines `Baselines` and the wavelenght `λ`. 
The first dimension of the returned uv plane is 2 (u and v)
"""
uvplane(Baselines, λ) = uvplane_(eltype(Baselines),Baselines,λ)
uvplane(Baselines, λ::Number)=uvplane(Baselines, [λ])

function uvplane(BaselinesU,BaselinesV, λ)
	Baselines = vcat(reshape(BaselinesU, 1,:),reshape(BaselinesV, 1,:))
	uvplane(Baselines, λ)
end

"""
    img = dirtybeam(uv,fov, pixsize)

Compute the dirtybeam image for measurements at uv coordinates `uv` spanning on the field of view `fov` with pixel size `pixsize`
"""
function dirtybeam(uv,fov, pixsize)
	N = round(Int,fov/pixsize)
	scaleduv = hcat([0; 0], reshape(uv* pixsize,2,:))  .|> NoUnits
	nfftplan  = plan_nfft(scaleduv,(N,N));
	δ = fill(1.0 /(N^2) + 0im,size(scaleduv,2))
	δ[1,1] = 1/(N^2) 
	beam = nfftplan' * δ
	return real.(beam)
end



"""
rx,ry,θ = beamellipse(uvp) 

Compute the parameters of the ellipse fitted on the central lobe of the dirty beam where:
- rx  is the semi-major axis 
- ry is the semi-minor axis
- θ is the principal angle
"""
function beamellipse(uvp::AbstractArray{T}) where {T}
	N = length(uvp[1,..])
	if T<:Quantity
		uvp = ustrip.(u"rad^-1",uvp)
	end
	S11 = 2 * sum(uvp[1,..].^2 ) / (2N+1)
	S22 = 2 * sum(uvp[2,..].^2 ) / (2N+1)
	S12 = 2 * sum(uvp[1,..] .* uvp[2,..] ) / (2N+1)
	#S = SHermitianCompact{2,Float64}([ S11 S12; S12 S22])
	S = [ S11 S12; S12 S22]
	vals, vecs = eigen(S)
	rx,ry = sqrt(2*log(2))/(2*π) ./ sqrt.(vals) .* u"rad" .|> mas
	θ  = atan(u"rad",vecs[1,1],vecs[2,1])
	return rx,ry, θ
end





function buildcovariance(rx,ry,θ)
	Vx,Vy = ((rx,ry) ./(2*sqrt(2*log(2)))).^2
	R =  [ 	cos(θ)  -sin(θ) ;
			sin(θ)  cos(θ) ]
	S =  [ 	Vx  0  ;
			0  	Vy ]
	return R'*S*R
end

function gaussian2D(tx,ty,W)
	r = tx.^2 * W[1,1] .- 2*W[1,2] *tx*ty' .+ (ty'.^2)* W[2,2]
	return 1/(2π) * sqrt(det(W)) .* exp.(- 0.5 .* r)
end


function __init__()
    @require Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80" begin
		#using Plots
		function plotimage(image, fov;  kw...)
			sz = size(image)
			Plots.heatmap(fftshift(fftfreq(sz[1])).*fov,fftshift(fftfreq(sz[1])).*fov,image , ticks=:native,aspect_ratio=:equal; kw...)
		end

		function plotuv(uv; kw...)
			Plots.plot(vcat(uv[1,:,:],-uv[1,:,:]),vcat(uv[2,:,:],-uv[2,:,:]), seriestype=:scatter,ticks=:native,aspect_ratio=:equal; kw...)
		end

		plotbeam!(p,rx,ry,θ; kw...) = Plots.plot!(p,t-> rx * cos(θ) * cos(t)  - ry * sin(θ) * sin(t), t-> rx * sin(θ) * cos(t)  + ry * cos(θ) * sin(t), 0, 2π, line=4,leg=false,ticks=:native; kw...)
			
    end
end
end # module BeamEstimation
