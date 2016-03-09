Subroutine FWald_quantitative(X, y, n, d, test_stat)
  Implicit None
  Integer, Intent(In) :: n,d
  Double Precision, Dimension(n), Intent(In) :: y
  Double Precision, Dimension(n,d), Intent(In) :: X
  Double Precision, Dimension(d) :: test_stat
  Double Precision :: YnormSq, kappa, r
  Double Precision, Dimension(n) :: Y_Ybar, X_Xbar
  Integer :: j

  Y_Ybar = y - Sum(y)/n
  YnormSq = Sum( Y_Ybar**2 )
  kappa = sqrt( DBLE(n-2) )
  Do j = 1,d
    X_Xbar = X(:,j) - Sum( X(:,j) )/n
    r = Dot_product(X_Xbar,Y_Ybar) / Sqrt( Sum( X_Xbar**2 ) * YnormSq )
    test_stat(j) = r * kappa/Sqrt(1-r**2)
  End Do
End Subroutine FWald_quantitative

Subroutine Fscore_qualitative(X, y, n, d, test_stat)
  Implicit None
  Integer, Intent(In) :: n,d
  Double Precision, Dimension(n), Intent(In) :: y
  Double Precision, Dimension(n,d), Intent(In) :: X
  Double Precision, Dimension(d) :: test_stat
  Double Precision :: U, V, Ymean
  Double Precision, Dimension(n) :: Y_Ybar
  Integer :: j

  test_stat = 0;   ! Default to 0
  Ymean = Sum(y)/n
  Y_Ybar = y - Ymean
  Do j = 1,d
    U = Dot_product( X(:,j), Y_Ybar )
    V = Ymean * (1-Ymean) * Sum( ( X(:,j) - Sum(X(:,j))/n )**2 )
    test_stat(j) = U**2 / V
  End Do
End Subroutine Fscore_qualitative

Subroutine Fscore_quantitative(X, y, n, d, test_stat)
  Implicit None
  Integer, Intent(In) :: n,d
  Double Precision, Dimension(n), Intent(In) :: y
  Double Precision, Dimension(n,d), Intent(In) :: X
  Double Precision, Dimension(d) :: test_stat
  Double Precision :: U, V, sigmaSq
  Double Precision, Dimension(n) :: Y_Ybar
  Integer :: j

  test_stat = 0;   ! Default to 0
  Y_Ybar = y - Sum(y)/n
  Do j = 1,d
    U = Dot_product( X(:,j), Y_Ybar )
    sigmaSq = Sum( (Y_Ybar)**2 ) / (n-1)
    V = sigmaSq*Sum( ( X(:,j) - Sum(X(:,j))/n )**2 )
    test_stat(j) = U**2 / V
  End Do
End Subroutine Fscore_quantitative

