import Whirl: Fields
using Fields

@testset "Point-Field Routines" begin

  @testset "Point creation" begin
    @test_throws AssertionError Points([1,2,3],[1,2])

    f = Points(
  end



end
