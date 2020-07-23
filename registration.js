
const store = new Vuex.Store({
    state:{
        user: {},
        expiration: 0
      },      
      mutations: {
        add_user(state, user) {
          state.user = user
          },
          set_expiration(state, expiration) {
              state.expiration = expiration
          },
          reset_state(state) {
            state.user = {}
            state.expiration = 0
          }
    },
    getters: {
        getExpiration() {
            return state.expiration;
        }
    }
      
})
  
Vue.component('Registration', {
    template: `
<v-dialog persistent v-model="notregistered" width="1500">
<v-card>
<v-form
id="registration-form"
v-model="valid"
@submit.prevent="processForm"
>
<v-container v-if="!this.$store.state.user.userId">
  <v-row justify="center">
    <v-col cols="10" md="6" lg="5" justify="center">
      <h2>Please Register as a FAUST User</h2>
      <v-text-field
        v-model="firstname"
        :rules="nameRules"
        label="First name"
        required
      ></v-text-field>
      <v-text-field
        v-model="lastname"
        :rules="nameRules"
        label="Last name"
        required
      ></v-text-field>

      <v-text-field
        v-model="email"
        :rules="emailRules"
        label="E-mail"
        required
      ></v-text-field>
      <v-text-field
        v-model="institution"
        :rules="institutionRules"
        label="Institution"
        required
      ></v-text-field>

      <v-btn
        :disabled="!valid"
        type="submit"
        form="registration-form"
        @click="hideModal"
        >Submit</v-btn
      >
    </v-col>
  </v-row>
</v-container>
</v-form>
</v-card >
</v-dialog>
`,
    data: () => ({
        notregistered: true,
        valid: false,
        firstname: '',
        lastname: '',
        nameRules: [
            (v) => !!v || 'Name is required',
            (v) => v.length <= 10 || 'Name must be less than 10 characters',
        ],
        email: '',
        emailRules: [
            (v) => !!v || 'E-mail is required',
            (v) => /.+@.+/.test(v) || 'E-mail must be valid',
        ],
        institution: '',
        institutionRules: [(v) => !!v || 'Institution is required'],
    }),
    mounted() {
        this.showModal()
        if (localStorage.user && localStorage.expiration)
        {
            this.hideModal()
            document.getElementById("navbuttons").classList.toggle("hidden")
            var user = JSON.parse(localStorage.user);
            this.userId = user.userId;
            this.institution = user.institution;
            this.email = user.email;
            this.firstname = user.name.split(" ")[0];
            this.lastname = user.name.split(" ")[1];
            this.$store.commit('add_user', user)
            this.$store.commit('set_expiration', localStorage.expiration)
            if (this.expired_min_old > 30)
            {
                localStorage.clear()
                this.userId = ''
                this.institution = ''
                this.email = ''
                this.firstname = ''
                this.lastname = ''
                this.$store.commit('reset_state')
                this.showModal()
                document.getElementById("navbuttons").classList.toggle("hidden")
            }
        }
    },
    computed: {
        expired_min_old() {
            return ((Date.now() - this.$store.state.expiration) / 1000 / 60 / 60 / 24)
        }
    },
    methods: {
        showModal() {
            this.notregistered = true;
        },
        hideModal() {
            this.notregistered = false;
        },
        processForm() {
            axios
                .post(
                    'https://fl50xsrcnf.execute-api.us-west-2.amazonaws.com/dev/users',
                    {
                        userId: Date.now().toString(),
                        name: this.firstname + ' ' + this.lastname,
                        email: this.email,
                        institution: this.institution,
                    }
                )
                .then((res) => {
                    this.$store.commit('add_user', res.data);
                    localStorage.user = JSON.stringify(res.data);
                    var expiration = Date.now()
                    this.$store.commit('set_expiration', expiration)
                    localStorage.expiration = expiration
                }).then(() => {
                    this.hideModal()
                    document.getElementById("navbuttons").classList.toggle("hidden")
                })
                .catch((err) => {
                    console.log(err)
                    this.showModal()
                });
    },
  },
});
const vueApp = new Vue({
  el: '#app',
  data: {
    display: 'redbox',
    },
    store,
    vuetify: new Vuetify()
  
});
