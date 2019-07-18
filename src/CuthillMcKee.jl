# Taken from https://github.com/rleegates/CuthillMcKee.jl/blob/master/src/CuthillMcKee.jl

# All credit goes to rleegates

using SparseArrays


"""
	 rcmpermute(A)
Computes and applies the reverse Cuthill-McKee permutation to the structurally-symmetric sparse matrix `A`. No checks are performed to ensure that the matrix is indeed structurally symmetric.
"""
function rcmpermute(A::SparseMatrixCSC)
	adj = coladjacency(A)
	sdeg = sortperm(coldegrees(A))
	p = symrcm(adj, sdeg)
	return A[p,p]
end

"""
	 symrcm(A,[ rev, inv, warnunconnected])
Computes the Cuthill-McKee permutation of the structurally-symmetric sparse matrix `A`. The default is to compute the reverse permutation, i.e. `rev` defaults to `true`. Optionally, when `inv` is set to `true`, the permutation is inverted. Per default the warning displaying the number of unconnected regions in the graph is disabled, to enable, set `warnunconnected` to `true`. No checks are performed to ensure that the matrix is indeed structurally symmetric.
"""
function symrcm(A::SparseMatrixCSC, args...)
	adj = coladjacency(A)
	sdeg = sortperm(coldegrees(A))
	return symrcm(adj, sdeg, args...)
end

"""
	 symrcm(adjac, sdegr,[ rev, inv, warnunconnected])
Computes the Cuthill-McKee permutation of the graph given as the adjacency list `adjac`. The adjacency list is required to be of type `Vector{Vector{Int}}` such that `adjac[i]` is a vector of the adjacent vertex numbers to vertex `i` in order of ascending vertex degree. The additional required argument `sdegr` is an array of vertex numbers in order of ascending vertex degree. The default is to compute the reverse permutation, i.e. `rev` defaults to `true`. Optionally, when `inv` is set to `true`, the permutation is inverted. Per default the warning displaying the number of unconnected regions in the graph is disabled, to enable, set `warnunconnected` to `true`.
"""
function symrcm(adjac::Array{Array{Int64,1}}, sdegr::Array{Int64,1}, rev=true, invert=false, warnunconnected=false)
	# Based on the description of the RCM algorithm by Ciprian Zavoianu.
	# See: http://ciprian-zavoianu.blogspot.com/2009/01/project-bandwidth-reduction.html
	nverts = length(adjac)
	checkbounds(sdegr,nverts)
	sd = copy(sdegr)
	R = Array{Int}(undef, 0)
	Q = Array{Int}(undef, 0)
	P = 0
	C = 0
	unconnected = 0
	inserted = falses(nverts)
	while !isempty(sd)
		found = false
		while !isempty(sd)
			P = popfirst!(sd)
			if !inserted[P]
				push!(R, P)
				inserted[P] = true
				found = true
				break
			end
		end
		if found
			empty!(Q)
			append!(Q, adjac[P])
			while !isempty(Q)
				while !isempty(Q)
					C = popfirst!(Q)
					if !inserted[C]
						push!(R, C)
						inserted[C] = true
						append!(Q, adjac[C])
						break
					end
				end
			end
		end
		if length(R) == nverts
			break
		else
			unconnected += 1
		end
	end
	if unconnected > 0 && warnunconnected
		warn("There are $(unconnected) unconnected regions.")
	end
	if rev
		reverse!(R)
	end
	if invert
		L = zeros(Int64,nverts)
		for i = 1:nverts
			L[R[i]] = i
		end
		return L
	else
		return R
	end
end

function coldegrees(A::SparseMatrixCSC)
	# assumes A is structurally symmetric
	cptr = A.colptr
	ncols = length(cptr)-1
	degr = zeros(Int,ncols)
	@inbounds for j = 1:ncols
		degr[j] = coldeg(A, j)
	end
	return degr
end

function coldeg(A::SparseMatrixCSC, j::Int)
	# assumes A is structurally symmetric
	cptr = A.colptr
	return cptr[j+1]-cptr[j]
end

function coladjacency(A::SparseMatrixCSC)
	# assumes A is structurally symmetric
	cptr = A.colptr
	rval = A.rowval
	ncols = length(cptr)-1
	adjac = Vector{Vector{Int}}(undef, ncols)
	sbyf = let A = A
		j->coldeg(A, j)
	end
	for j = 1:ncols
		strt = cptr[j]
		jdeg = coldeg(A, j)
		jadj = Vector{Int}(undef, jdeg)
		for i = 1:jdeg
			jadj[i] = rval[strt+i-1]
		end
		sort!(jadj, by=sbyf)
		adjac[j] = jadj
	end
	return adjac
end
