struct VecArray{s_, T, dim}
    data::Array{T, dim}
end

Base.eltype(::Type{VecArray{s, T, dim}}) where {s, T, dim} = T

@inline function Base.getindex(v::VecArray{s, T, dim}, k...) where {s, T, dim}
    #    Vec{s, T}(NTuple{s, T}(@inbounds v.data[is, k...] for is in 1:s))
    return Vec{s, T}(ntuple(is -> @inbounds(v.data[is, k...]), s))
end

@inline function Base.setindex!(
        v::VecArray{s, T, dim}, vi::Vec{s, T}, k...
    ) where {s, T, dim}
    @inbounds for is in 1:s
        v.data[is, k...] = vi[is]
    end
    return nothing
end

@inline function Base.setindex!(v::VecArray{s, T, dim}, vk::T2, k...) where {s, T, T2, dim}
    vk_ = convert(T, vk)
    @inbounds for is in 1:s
        v.data[is, k...] = vk_
    end
    return nothing
end

## getindex_ and  setindex_! implementations

@inline function getindex(v::VecArray{s, T, dim}, k::Int64) where {s, T, dim}
    j = s * (k - 1)
    #    Vec{s, T}(NTuple{s, T}(@inbounds v.data[is + j] for is in 1:s))
    return Vec{s, T}(ntuple(is -> @inbounds(v.data[is + j]), s))
end

@inline function setindex!(
        v::VecArray{s, T, dim}, vk::Vec{s, T}, k::Int64
    ) where {s, T, dim}
    j = s * (k - 1)
    @inbounds for is in 1:s
        v.data[is + j] = vk[is]
    end
    return nothing
end

@inline function setindex!(v::VecArray{s, T, dim}, vk::T2, k::Int64) where {s, T, T2, dim}
    vk_ = convert(T, vk)
    j = s * (k - 1)
    @inbounds for is in 1:s
        v.data[is + j] = vk_
    end
    return nothing
end
