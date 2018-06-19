#include <dune/geometry/quadraturerules.hh>

template< class Surface >
double calcRadius( const Surface &surface )
{
  typedef typename Surface::template Codim< 0 >::Entity Entity;

  double R = 0.;
  double vol = 0.;
  for( const Entity &entity : elements( surface ) )
  {
    const auto& rule = Dune::QuadratureRules<double, 2>::rule(entity.type(), 4);
    for (const auto& p : rule)
    {
      const auto geo = entity.geometry();
      const double weight = geo.volume() * p.weight();
      R   += geo.toGlobal(p.position()).twoNorm() * weight;
      vol += weight;
    }
    return R/vol;
}

/*
# compute an averaged radius of the surface
def calcRadius(surface):
    # compute R = int_x |x| / int_x 1
    R   = 0
    vol = 0
    for e in surface.elements:
        rule = geometry.quadratureRule(e.type, 4)
        for p in rule:
            geo = e.geometry
            weight = geo.volume * p.weight
            R   += geo.toGlobal(p.position).two_norm * weight
            vol += weight
    return R/vol
*/
