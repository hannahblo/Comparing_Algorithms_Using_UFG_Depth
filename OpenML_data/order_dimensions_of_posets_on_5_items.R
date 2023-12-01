# This file computes the order dimensions for all posets on 5 items

context_OpenML <- ddandrda:::compute_context_all_p_orders(n_items=5)
all_posets_5 <- ddandrda:::compute_all_partial_orders(n_items=5,list=TRUE,complemented=FALSE)
order_dimensions_5 <- rep(0,4231)

for(k in (1:4231)){
    order_dimensions_5[k] <- oofos:::compute_order_dimension(all_posets_5[[k]],context_OpenML)
	print(order_dimensions_5[k])
}


max(order_dimensions_5)
# [1] 2
# The maximal order dimension of a poset on 5 items is 2
