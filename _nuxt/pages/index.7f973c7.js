(window.webpackJsonp=window.webpackJsonp||[]).push([[3],{307:function(t,e,n){"use strict";n.r(e);n(16),n(17),n(4);var r=n(272),l=n.n(r),o={data:function(){return{valid:!1,firstname:"",lastname:"",nameRules:[function(t){return!!t||"Name is required"},function(t){return t.length<=10||"Name must be less than 10 characters"}],email:"",emailRules:[function(t){return!!t||"E-mail is required"},function(t){return/.+@.+/.test(t)||"E-mail must be valid"}],institution:"",institutionRules:[function(t){return!!t||"Institution is required"}]}},methods:{processForm:function(){var t=this;l.a.post("https://fl50xsrcnf.execute-api.us-west-2.amazonaws.com/dev/users",{userId:Date.now().toString(),name:this.firstname+" "+this.lastname,email:this.email,institution:this.institution}).then((function(t){return console.log(t),t})).then((function(e){t.$store.commit("user/add_user",e.data)})).then(this.$router.push("/success")).catch((function(t){return console.log(t)}))}}},c=n(64),m=n(93),d=n.n(m),f=n(254),v=n(270),h=n(256),x=n(304),w=n(271),_=n(305),component=Object(c.a)(o,(function(){var t=this,e=t.$createElement,n=t._self._c||e;return n("v-container",[n("h2",[t._v("Register to Download FAUST")]),t._v(" "),n("v-form",{attrs:{id:"registration-form"},on:{submit:function(e){return e.preventDefault(),t.processForm(e)}},model:{value:t.valid,callback:function(e){t.valid=e},expression:"valid"}},[n("v-container",[n("v-row",[n("v-col",{attrs:{cols:"12",md:"3"}},[n("v-text-field",{attrs:{rules:t.nameRules,label:"First name",required:""},model:{value:t.firstname,callback:function(e){t.firstname=e},expression:"firstname"}})],1),t._v(" "),n("v-col",{attrs:{cols:"12",md:"3"}},[n("v-text-field",{attrs:{rules:t.nameRules,label:"Last name",required:""},model:{value:t.lastname,callback:function(e){t.lastname=e},expression:"lastname"}})],1),t._v(" "),n("v-col",{attrs:{cols:"12",md:"3"}},[n("v-text-field",{attrs:{rules:t.emailRules,label:"E-mail",required:""},model:{value:t.email,callback:function(e){t.email=e},expression:"email"}})],1),t._v(" "),n("v-col",{attrs:{cols:"12",md:"3"}},[n("v-text-field",{attrs:{rules:t.institutionRules,label:"Institution",required:""},model:{value:t.institution,callback:function(e){t.institution=e},expression:"institution"}})],1),t._v(" "),n("v-btn",{staticClass:"color: primary",attrs:{disabled:!t.valid,type:"submit",form:"registration-form"}},[t._v("Submit")])],1)],1)],1)],1)}),[],!1,null,null,null);e.default=component.exports;d()(component,{VBtn:f.a,VCol:v.a,VContainer:h.a,VForm:x.a,VRow:w.a,VTextField:_.a})}}]);