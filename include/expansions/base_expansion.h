#ifndef BASE_EXPANSION_H
#define BASE_EXPANSION_H

class base_expansion {
public:
    virtual ~base_expansion() = 0;

    [[nodiscard]] virtual bool validate() const = 0;

    virtual void apply() const = 0;
};

#endif //BASE_EXPANSION_H
