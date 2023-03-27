#ifndef _levicivita_h_
#define _levicivita_h_

double leviCivita(const double& i,
                  const double& j,
                  const double& k)
{
    return 0.5*(i-j)*(j-k)*(k-i);
}

#endif  //_levicivita_h_
