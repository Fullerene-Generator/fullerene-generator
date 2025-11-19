#include<expansions/f_expansion.h>

bool f_expansion::validate() const {
    for (const auto& u : v_->neighbors()) {
        if (u.lock()->type() == node_type::NODE_6) {
            return false;
        }
    }

    return true;
}

void f_expansion::apply() const {

}