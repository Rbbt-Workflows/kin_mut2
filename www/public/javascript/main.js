function show_hide(){
  var target = $("#" + $(this).attr('ref'));
  target.toggle();
  return false;
}

function tabs(elem){
 var tabs = elem.find('li > a')
 tabs.click(select_tab)
 
 var active = elem.find('li > a.active').first();

 if (active == null)
  tabs.first().click();
 else
  active.click();

}

function select_tab(){
  var target = $(this).attr('target');
  var tab_list = $(this).parent('li').parent('ul');

  if ($(this).attr('class') == 'disabled') return(false) 

  var other_tabs = tab_list.find('li > a');

  other_tabs.each(function(key,value){
    var other_target = $(value).attr('target');
    $(value).removeClass('active')
    $(other_target).hide()
  })

  $(target).show()
  $(this).addClass('active')


  return(false);
}


