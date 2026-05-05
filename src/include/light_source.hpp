#ifndef LIGHT_SOURCE_HPP
#define LIGHT_SOURCE_HPP

class LightSource {
  public:
    int source_type;
    double limb[4];
    double normalization;

    LightSource();
    LightSource(int type, double limb0, double limb1, double limb2, double limb3);

    double get_brightness(double x, double y);
};

#endif
