using BenchmarkTools

struct non_typed
    a
    b
    c
end

struct typed
    a::Int64
    b::Int64
    c::Int64
end

struct parametric{T}
    a::T
    b::T
    c::T
end

struct inception{T,S}
    a::T
    b::T
    c::S
end

function addstuff(a)
    return a.a + a.b
end

function addstuff_inception(a)
    return a.a.a + a.b.a
end


a = non_typed(1,2,3)
b = typed(1,2,3)
c = parametric(1,2,3)


@code_warntype addstuff(a)
@code_warntype addstuff(b)
@code_warntype addstuff(c)

@btime addstuff(a)
@btime addstuff(b)
@btime addstuff(c)


struct tts{T,S}
    a::T
    b::T
    c::S
end

struct asdf{T,S}
    a::T
    b::T
    c::S
    function asdf(s,t,u)
        a = s
        b = t
        c = u
        return new{typeof(a),typeof(b)}(a,b,c)
    end
end
