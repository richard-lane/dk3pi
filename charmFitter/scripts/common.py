
def reset_params(Fitter, decay_params):
    Fitter.setParameter("x", decay_params[0])
    Fitter.setParameter("y", decay_params[1])
    Fitter.setParameter("r", decay_params[2])
    Fitter.setParameter("z_im", decay_params[3])
    Fitter.setParameter("z_re", decay_params[4])
    Fitter.setParameter("width", decay_params[5])

