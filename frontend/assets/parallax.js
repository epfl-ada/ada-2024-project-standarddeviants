var $ = jQuery.noConflict();
$(window).scroll(function(){
  var scroll = $(this).scrollTop();
  $(".parallax").css({
    "background-position-x": "center",
    "background-position-y": scroll/2 + "px"
  });
});
