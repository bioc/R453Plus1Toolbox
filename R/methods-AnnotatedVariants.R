setMethod("annotatedVariants", signature(object="AnnotatedVariants"),
    function(object) {
        return(object@annotatedVariants)
    }
)

setReplaceMethod("annotatedVariants",
    signature=signature(object="AnnotatedVariants", value="list"),
    function(object, value) {
        object@annotatedVariants = value
        return(object)
    }
)

setMethod("names", signature(x="AnnotatedVariants"),
    function(x) {
        return(names(x@annotatedVariants))
    }
)

setReplaceMethod("names",
    signature=signature(x="AnnotatedVariants", value="character"),
    function(x, value) {
        names(x@annotatedVariants) = value
        return(x)
    }
)

setMethod("show", signature(object="AnnotatedVariants"),
    function(object) {
        cat(paste("An object of class AnnotatedVariants with",
            length(annotatedVariants(object)), "variants.\n"))
    }
)
