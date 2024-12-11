var $ = jQuery.noConflict()
$(window).scroll(function(){
  var scroll = $(this).scrollTop()
  $(".parallax").css({"background-position":"0px "+scroll/2+"px"})
})
