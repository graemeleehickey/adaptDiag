# binom_sample_size errors when p0 >= p1

    Code
      binom_sample_size(p0 = 0.9, p1 = 0.8)
    Condition
      Error in `binom_sample_size()`:
      ! p0 must be less than p1

---

    Code
      binom_sample_size(p0 = 0.8, p1 = 0.8)
    Condition
      Error in `binom_sample_size()`:
      ! p0 must be less than p1

