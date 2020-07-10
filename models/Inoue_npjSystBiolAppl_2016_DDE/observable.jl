const observables = [
    "Nuclear_NFkB"
]

function observables_index(observable_name::String)::Int

    return findfirst(isequal(observable_name),observables)
end