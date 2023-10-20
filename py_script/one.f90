module one
use two

contains

function onefunc(x) result(y)
    real(8) :: x,y
    y = twofunc(x)
    y = twofunc(y)
end function

end module