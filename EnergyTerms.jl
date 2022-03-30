# ENERGY TERMS
# to be used with the module QUANTUM DYNAMICS
# by Jonas Boym Flaten


function within(groups::Vector{Int64};bodies::Int64=2)
    # returns a vector of tuples of the given particle group indices with the given number of bodies.
    return [fill(group,bodies) for group in groups]
end
within(groups...;bodies=2) = within([groups...];bodies=bodies)

function between(groups::Vector{Int64};bodies::Int64=2)
    # returns a vector of all combinations from the given particle groups with the given number of bodies.
    @assert (bodies ≤ length(groups)) "The number of involved bodies is larger than the number of given particle groups; "*
        "input at least as many indices in 'groups' as 'bodies'."
    if (bodies == 1)
        return [[group] for group in groups]
    else
        combinations = Vector{Int64}[]
        for i in 1:(length(groups)-(bodies-1))
            for combination in between(groups[(i+1):end];bodies=(bodies-1))
                push!(combinations,[groups[i],combination...])
            end
        end
        return combinations
    end
end
between(groups...;bodies=2) = between([groups...];bodies=bodies)
bethreen(groups) = between(groups;bodies=3)
bethreen(groups...) = bethreen([groups...])


abstract type EnergyTerm end
    # comprises various possible energy terms of the system Hamiltonian.

@enum ParticleProperty m q s r p l # encodes the particle properties involved in energy terms.

struct CustomEnergyTerm <: EnergyTerm
    # represents arbitrary energy terms for the given particle groups, specified by involved properties and form.
    groups::Vector{Vector{Int64}} # are the combinations of particle groups involved with the energy term.
    properties::Vector{Vector{ParticleProperty}} # are the particle properties involved in the energy term, sorted by body.
    form::Function # is the mathematical form of the energy term, declared in the same order as the particle properties above.
    function CustomEnergyTerm(groups,properties,form)
        bodies = length(properties) # is the number of bodies involved in the energy term.
        for combination in groups
            @assert (length(combination) == bodies) "The number of involved bodies does not match the given particle group combinations; "*
                "input tuples in 'groups' whose sizes match the number of tuples in 'properties'."
        end
        form_arguments = first(methods(form)).nargs-1 # is the number of arguments in the mathematical form of the energy term.
        @assert (form_arguments == sum([length(body) for body in properties])) "The declared number of arguments "*
            "does not match the declared number of properties; input as many elements in 'properties' as there are arguments in 'form'."
        return new(groups,properties,form)
    end
end
CustomEnergyTerm(groups;properties,form) = CustomEnergyTerm(groups,properties,form)
function CustomEnergyTerm(within_groups::Vector{Int64},between_groups::Vector{Int64},properties,form)
    # returns a CustomEnergyTerm which applies within the particle groups given in 'within_groups'
        # and between the particle groups given in 'between_groups'.
    bodies = length(properties) # is the number of bodies involved in the energy term.
    groups = [within(within_groups)...,between(between_groups;bodies=bodies)...]
    return CustomEnergyTerm(groups,properties,form)
end
CustomEnergyTerm(within_groups,between_groups;properties,form) = CustomEnergyTerm(within_groups,between_groups,properties,form)
CustomEnergyTerms(properties,form) = CustomEnergyTerm([fill(0,length(properties))],properties,form)
    # returns a CustomEnergyTerm applying to all particles in the system. (groups = [0] encodes that the term applies to all particles.)
CustomEnergyTerms(;properties,form) = CustomEnergyTerms(properties,form)

struct KineticTerm <: EnergyTerm
    # represents one-body kinergy terms for the given particle groups, of the form p²/2m.
    groups::Vector{Int64} # are the particle groups involved with a kinetic term.
end
KineticTerm(group) = KineticTerm([group])
KineticTerm(groups...) = KineticTerm([groups...])
KineticTerms() = KineticTerm([0])
    # returns a KineticTerm applying to all particles in the system. (groups = [0] encodes that the term applies to all particles.)

struct PotentialTerm <: EnergyTerm
    # represents arbitrary one-body potergy terms for the given particle groups, specified by involved properties and form.
    groups::Vector{Int64} # are the particle groups involved with the potential term.
    properties::Vector{ParticleProperty} # are the particle properties involved in the potential term.
    form::Function # is the mathematical form of the potential term, declared in the same order as the particle properties above.
    function PotentialTerm(groups,properties,form)
        valid_properties = (m,q,s,r)
        for property in properties
            @assert (property ∈ valid_properties) "The particle property $property is not allowed in a potential energy term; "*
                "input only valid properties $valid_properties."
        end
        form_arguments = first(methods(form)).nargs-1 # is the number of arguments in the mathematical form of the potential term.
        @assert (form_arguments == length(properties)) "The declared number of arguments "*
            "does not match the declared number of properties; input as many elements in 'properties' as there are arguments in 'form'."
        return new(groups,properties,form)
    end
end
PotentialTerm(groups;properties,form) = PotentialTerm(groups,properties,form)
PotentialTerms(properties,form) = PotentialTerm([0],properties,form)
    # returns a PotentialTerm applying to all particles in the system. (groups = [0] encodes that the term applies to all particles.)
PotentialTerms(;properties,form) = PotentialTerms(properties,form)

struct HarmonicTerm <: EnergyTerm
    # represents one-body harmonic potergy terms for the given particle groups,
        # of the form m(ω²/2)(r-r₀)² with a strength ω²/2 and a centre r₀.
    groups::Vector{Int64} # are the particle groups involved with the harmonic term.
    strength::Float64 # is the strength ω²/2 of the harmonic term.
    centre::Vector{Float64} # is the centre r₀ of the harmonic term.
    function HarmonicTerm(groups,strength=1.0,centre=[0.,0.,0.])
        @assert (strength > 0) "The given harmonic strength $strength is invalid; input a positive number as 'strength'."
        return new(groups,strength,centre)
    end
end
HarmonicTerm(groups;strength,centre=[0.,0.,0.]) = HarmonicTerm(groups,strength,centre)
HarmonicTerms(strength=1.0,centre=[0.,0.,0.]) = HarmonicTerm([0],strength,centre)
    # returns a HarmonicTerm applying to all particles in the system. (groups = [0] encodes that the term applies to all particles.)
HarmonicTerms(;strength,centre=[0.,0.,0.]) = HarmonicTerms(strength,centre)

struct Newton1Term <: EnergyTerm
    # represents one-body Newton potergy terms for the given particle groups,
        # of the form Am/√((r-r₀)²+a²) with a strength A, a centre r₀ and a shielding a.
    groups::Vector{Int64} # are the particle groups involved with the Newton term.
    strength::Float64 # is the strength A of the Newton term.
    centre::Vector{Float64} # is the centre r₀ of the Newton term.
    shielding::Float64 # is the shielding a of the Newton term.
    function Newton1Term(groups,strength=1.0,centre=[0.,0.,0.],shielding=0.1)
        @assert (strength > 0) "The given Newton strength $strength is invalid; input a positive number as 'strength'."
        return new(groups,strength,centre,shielding)
    end
end
Newton1Term(groups;strength,centre=[0.,0.,0.],shielding=0.1) = Newton1Term(groups,strength,centre,shielding)
Newton1Terms(strength=1.0,centre=[0.,0.,0.],shielding=0.1) = Newton1Term([0],strength,centre,shielding)
    # returns a Newton1Term applying to all particles in the system. (groups = [0] encodes that the term applies to all particles.)
Newton1Terms(;strength,centre=[0.,0.,0.],shielding=0.1) = Newton1Terms(strength,centre,shielding)

struct Newton2Term <: EnergyTerm
    # represents two-body Newton interaction terms for the given particle groups,
        # of the form Am₁m₂/√((r₁-r₂)²+a²) with a strength A and a shielding a.
    groups::Vector{NTuple{2,Int64}} # are the combinations of particle groups involved with Newton interaction.
    strength::Float64 # is the strength A of the Newton term.
    shielding::Float64 # is the shielding a of the Newton term.
    function Newton2Term(groups,strength,shielding)
        @assert (strength > 0) "The given Newton interaction strength $strength is invalid; input a positive number as 'strength'."
        return new(groups,strength,shielding)
    end
end
Newton2Term(groups;strength=1.0,shielding=0.1) = Newton2Term(groups,strength,shielding)
function Newton2Term(within_groups::Vector{Int64},between_groups::Vector{Int64},strength::Float64=1.0,shielding::Float64=0.1)
    # returns a Newton2Term which applies within the particle groups given in 'within_groups'
        # and between the particle groups given in 'between_groups'.
    groups = [within(within_groups)...,between(between_groups;bodies=2)...]
    return Newton2Term(groups,strength,shielding)
end
Newton2Term(within_groups,between_groups;strength,shielding=0.1) = Newton2Term(within_groups,between_groups,strength,shielding)
Newton2Terms(strength=1.0,shielding=0.1) = Newton2Term([(0,0)],strength,shielding)
    # returns a Newton2Term applying to all particles in the system. (groups = [(0,0)] encodes that the term applies to all particles.)
Newton2Terms(;strength,shielding=0.1) = Newton2Terms(strength,shielding)

struct Coulomb1Term <: EnergyTerm
    # represents one-body Coulomb potergy terms for the given particle groups,
        # of the form Aq/√((r-r₀)²+a²) with a strength A, a centre r₀ and a shielding a.
    groups::Vector{Int64} # are the particle groups involved with the Coulomb term.
    strength::Float64 # is the strength A of the Coulomb term.
    centre::Vector{Float64} # is the centre r₀ of the Coulomb term.
    shielding::Float64 # is the shielding a of the Coulomb term.
    function Coulomb1Term(groups,strength=1.0,centre=[0.,0.,0.],shielding=0.1)
        @assert (strength > 0) "The given Coulomb strength $strength is invalid; input a positive number as 'strength'."
        return new(groups,strength,centre,shielding)
    end
end
Coulomb1Term(groups;strength,centre=[0.,0.,0.],shielding=0.1) = Coulomb1Term(groups,strength,centre,shielding)
Coulomb1Terms(strength=1.0,centre=[0.,0.,0.],shielding=0.1) = Coulomb1Term([0],strength,centre,shielding)
    # returns a Coulomb1Term applying to all particles in the system. (groups = [0] encodes that the term applies to all particles.)
Coulomb1Terms(;strength,centre=[0.,0.,0.],shielding=0.1) = Coulomb1Terms(strength,centre,shielding)

struct Coulomb2Term <: EnergyTerm
    # represents two-body Coulomb interaction terms for the given particle groups,
        # of the form Aq₁q₂/√((r₁-r₂)²+a²) with a strength A and a shielding a.
    groups::Vector{NTuple{2,Int64}} # are the combinations of particle groups involved with Coulomb interaction.
    strength::Float64 # is the strength A of the Coulomb term.
    shielding::Float64 # is the shielding a of the Coulomb term.
    function Coulomb2Term(groups,strength=1.0,shielding=0.1)
        @assert (strength > 0) "The given Coulomb interaction strength $strength is invalid; input a positive number as 'strength'."
        return new(groups,strength,shielding)
    end
end
Coulomb2Term(groups;strength,shielding=0.1) = Coulomb2Term(groups,strength,shielding)
function Coulomb2Term(within_groups::Vector{Int64},between_groups::Vector{Int64},strength::Float64=1.0,shielding::Float64=0.1)
    # returns a Coulomb2Term which applies within the particle groups given in 'within_groups'
        # and between the particle groups given in 'between_groups'.
    groups = [within(within_groups)...,between(between_groups;bodies=2)...]
    return Coulomb2Term(groups,strength,shielding)
end
Coulomb2Term(within_groups,between_groups;strength,shielding=0.1) = Coulomb2Term(within_groups,between_groups,strength,shielding)
Coulomb2Terms(strength=1.0,shielding=0.1) = Coulomb2Term([(0,0)],strength,shielding)
    # returns a Coulomb2Term applying to all particles in the system. (groups = [(0,0)] encodes that the term applies to all particles.)
Coulomb2Terms(;strength,shielding=0.1) = Coulomb2Terms(strength,shielding)

FreeHamiltonian() = [KineticTerms()]
HarmonicHamiltonian(strength=1.0,centre=[0.,0.,0.]) = [KineticTerms(),HarmonicTerms(strength,centre)]
HarmonicHamiltonian(;strength,centre=[0.,0.,0.]) = HarmonicHamiltonian(strength,centre)
NewtonHamiltonian(strength=1.0,shielding=0.1) = [KineticTerms(),Newton2Terms(strength,shielding)]
NewtonHamiltonian(;strength,shielding=0.1) = NewtonHamiltonian(strength,shielding)
CoulombHamiltonian(strength=1.0,shielding=0.1) = [KineticTerms(),Coulomb2Terms(strength,shielding)]
CoulombHamiltonian(strength,shielding=0.1) = CoulombHamiltonian(strength,shielding)

; # suppresses inclusion output.
